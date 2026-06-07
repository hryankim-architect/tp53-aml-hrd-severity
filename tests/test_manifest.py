"""Static guards on data/manifest.yaml + the checksum write-back logic.

These tests run offline. They assert the committed manifest is real (no inherited
scaffold stub, no zero/blank checksums) and that ``write_manifest_checksums``
fills hashes back correctly. The network paths (``fetch_manifest`` download,
``refresh_manifest`` GDC enumeration) are exercised by the real ``make data``
step, not unit-tested here.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

from tp53_hrd import gdc, pipeline

MANIFEST = Path(__file__).resolve().parents[1] / "data" / "manifest.yaml"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
STUB_SHA = "0" * 64


def _load() -> dict:
    return yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))


class TestManifestIsReal:
    def test_has_clinical_and_inputs(self):
        m = _load()
        assert "clinical" in m and "inputs" in m
        assert isinstance(m["inputs"], list) and len(m["inputs"]) > 0

    def test_no_scaffold_stub(self):
        text = MANIFEST.read_text(encoding="utf-8")
        # The inherited scaffold placeholder must be gone.
        assert "example.org" not in text
        assert STUB_SHA not in text

    def test_clinical_block_is_pinned(self):
        clin = _load()["clinical"]
        assert clin["endpoint"] == gdc.GDC_CASES_URL
        assert clin["path"] == "tcga-laml/clinical.json"
        assert HEX64.match(clin["sha256"]), "clinical sha256 must be a real 64-hex digest"

    def test_every_input_is_real_and_open_tier(self):
        inputs = _load()["inputs"]
        seen_paths: set[str] = set()
        for e in inputs:
            assert e["url"].startswith(gdc.GDC_DATA_URL + "/"), e["url"]
            assert e["path"].startswith("tcga-laml/mafs/")
            assert e["path"].endswith(".maf.gz")
            assert HEX64.match(e["sha256"]), f"input {e['path']} has a non-real sha256"
            assert e["sha256"] != STUB_SHA
            assert e["path"] not in seen_paths, f"duplicate path {e['path']}"
            seen_paths.add(e["path"])

    def test_url_and_path_uuid_agree(self):
        # path is <file_id>.maf.gz; url is .../data/<file_id> — they must match.
        for e in _load()["inputs"]:
            file_id = e["url"].rsplit("/", 1)[1]
            assert e["path"] == f"tcga-laml/mafs/{file_id}.maf.gz"


class TestWriteManifestChecksums:
    def test_fills_clinical_and_inputs_by_path(self, tmp_path):
        # A minimal manifest with blank checksums (post-refresh shape).
        src = (
            "clinical:\n"
            "  endpoint: https://api.gdc.cancer.gov/cases\n"
            "  path: tcga-laml/clinical.json\n"
            "  sha256:\n"
            "inputs:\n"
            "  - url: https://api.gdc.cancer.gov/data/aaa\n"
            "    path: tcga-laml/mafs/aaa.maf.gz\n"
            "    sha256:\n"
            "  - url: https://api.gdc.cancer.gov/data/bbb\n"
            "    path: tcga-laml/mafs/bbb.maf.gz\n"
            "    sha256:\n"
        )
        mpath = tmp_path / "manifest.yaml"
        mpath.write_text(src, encoding="utf-8")

        result = {
            "clinical": {"rel": "tcga-laml/clinical.json", "status": "downloaded", "sha256": "c" * 64},
            "inputs": [
                {"rel": "tcga-laml/mafs/aaa.maf.gz", "status": "downloaded", "sha256": "a" * 64},
                {"rel": "tcga-laml/mafs/bbb.maf.gz", "status": "downloaded", "sha256": "b" * 64},
            ],
        }
        n = pipeline.write_manifest_checksums(mpath, result)
        assert n == 3

        m = yaml.safe_load(mpath.read_text(encoding="utf-8"))
        assert m["clinical"]["sha256"] == "c" * 64
        assert m["inputs"][0]["sha256"] == "a" * 64
        assert m["inputs"][1]["sha256"] == "b" * 64

    def test_cached_results_write_nothing(self, tmp_path):
        mpath = tmp_path / "manifest.yaml"
        mpath.write_text("inputs:\n  - url: x\n    path: p\n    sha256:\n", encoding="utf-8")
        # 'cached' (not 'downloaded') results carry no new sha → nothing written.
        n = pipeline.write_manifest_checksums(
            mpath, {"inputs": [{"rel": "p", "status": "cached"}]}
        )
        assert n == 0


class TestBuildMafInputs:
    def test_shape_and_blank_sha(self):
        files = [
            {"file_id": "uuid-1", "file_name": "x.maf.gz", "file_size": 123, "md5sum": "m"},
        ]
        entries = gdc.build_maf_inputs(files)
        assert len(entries) == 1
        e = entries[0]
        assert e["url"] == gdc.GDC_DATA_URL + "/uuid-1"
        assert e["path"] == "tcga-laml/mafs/uuid-1.maf.gz"
        assert e["sha256"] == ""  # filled later by fetch --write-checksums
        assert e["size_bytes"] == 123
