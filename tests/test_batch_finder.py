#!/usr/bin/env python3
"""
Comprehensive tests for BatchFinder and family lookup functionality.

Tests the smart batch detection, protein searching, and configuration management.
"""

import os
import tempfile
from pathlib import Path

import pytest

from pyecod_mini.cli.config import BatchFinder, PyEcodMiniConfig


@pytest.mark.unit
class TestBatchFinderInitialization:
    """Test BatchFinder initialization and basic properties"""

    def test_batch_finder_initialization(self):
        """Test that BatchFinder can be initialized"""
        finder = BatchFinder("/tmp/test_batches")

        assert finder.base_dir == Path("/tmp/test_batches")
        assert isinstance(finder._batch_cache, dict)
        assert isinstance(finder.stable_batches, dict)

    def test_batch_finder_stable_batches(self):
        """Test that stable_batches registry exists"""
        finder = BatchFinder("/tmp/test")

        # Should have known test cases
        assert isinstance(finder.stable_batches, dict)
        assert "8ovp_A" in finder.stable_batches
        assert finder.stable_batches["8ovp_A"] == "ecod_batch_036_20250406_1424"

    def test_batch_finder_with_real_base_dir(self):
        """Test BatchFinder with real ECOD batch directory"""
        base_dir = "/data/ecod/pdb_updates/batches"

        if not Path(base_dir).exists():
            pytest.skip("ECOD batch directory not available")

        finder = BatchFinder(base_dir)
        assert finder.base_dir.exists()


@pytest.mark.unit
class TestBatchFinderAvailableBatches:
    """Test batch discovery functionality"""

    def test_get_available_batches_empty_dir(self):
        """Test _get_available_batches() with empty directory"""
        with tempfile.TemporaryDirectory() as tmpdir:
            finder = BatchFinder(tmpdir)
            batches = finder._get_available_batches()
            assert batches == []

    def test_get_available_batches_nonexistent_dir(self):
        """Test _get_available_batches() with non-existent directory"""
        finder = BatchFinder("/nonexistent/path")
        batches = finder._get_available_batches()
        assert batches == []

    def test_get_available_batches_with_mock_batches(self):
        """Test _get_available_batches() finds batch directories"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create mock batch directories
            Path(tmpdir, "ecod_batch_001_20250101_0000").mkdir()
            Path(tmpdir, "ecod_batch_002_20250201_0000").mkdir()
            Path(tmpdir, "alt_rep_batch_001_20250101_0000").mkdir()
            Path(tmpdir, "not_a_batch").mkdir()  # Should be ignored

            finder = BatchFinder(tmpdir)
            batches = finder._get_available_batches()

            # Should find ecod_batch and alt_rep_batch directories
            assert len(batches) >= 2
            assert any("ecod_batch_001" in b for b in batches)
            assert any("ecod_batch_002" in b for b in batches)
            assert any("alt_rep_batch_001" in b for b in batches)
            assert "not_a_batch" not in batches

            # Should be sorted
            assert batches == sorted(batches)

    @pytest.mark.integration
    def test_get_available_batches_real_data(self):
        """Test _get_available_batches() with real ECOD data"""
        base_dir = "/data/ecod/pdb_updates/batches"

        if not Path(base_dir).exists():
            pytest.skip("ECOD batch directory not available")

        finder = BatchFinder(base_dir)
        batches = finder._get_available_batches()

        # Should find batches
        assert len(batches) > 0

        # All should start with ecod_batch_ or alt_rep_batch_
        for batch in batches:
            assert batch.startswith("ecod_batch_") or batch.startswith("alt_rep_batch_")


@pytest.mark.unit
class TestBatchFinderProteinLookup:
    """Test protein existence checking"""

    def test_protein_exists_in_batch_false(self):
        """Test _protein_exists_in_batch() returns False when protein doesn't exist"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            (batch_dir / "domains").mkdir()

            finder = BatchFinder(tmpdir)
            exists = finder._protein_exists_in_batch("8abc_A", "ecod_batch_001_test")

            assert exists is False

    def test_protein_exists_in_batch_true(self):
        """Test _protein_exists_in_batch() returns True when protein exists"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()

            # Create mock domain summary file
            domain_file = domains_dir / "8ovp_A.develop291.domain_summary.xml"
            domain_file.write_text("<domain_summary/>")

            finder = BatchFinder(tmpdir)
            exists = finder._protein_exists_in_batch("8ovp_A", "ecod_batch_001_test")

            assert exists is True

    @pytest.mark.integration
    def test_protein_exists_in_batch_real_data(self):
        """Test _protein_exists_in_batch() with real ECOD data"""
        base_dir = "/data/ecod/pdb_updates/batches"

        if not Path(base_dir).exists():
            pytest.skip("ECOD batch directory not available")

        finder = BatchFinder(base_dir)

        # Check for stable test protein in its known batch
        exists = finder._protein_exists_in_batch("8ovp_A", "ecod_batch_036_20250406_1424")

        if Path(base_dir, "ecod_batch_036_20250406_1424").exists():
            # If batch exists, protein should be in it
            assert exists is True


@pytest.mark.unit
class TestBatchFinderGetProteins:
    """Test getting list of proteins in batch"""

    def test_get_proteins_in_batch_empty(self):
        """Test _get_proteins_in_batch() with empty batch"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            (batch_dir / "domains").mkdir()

            finder = BatchFinder(tmpdir)
            proteins = finder._get_proteins_in_batch("ecod_batch_001_test")

            assert proteins == []

    def test_get_proteins_in_batch_with_proteins(self):
        """Test _get_proteins_in_batch() finds proteins"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()

            # Create mock domain summary files
            (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")
            (domains_dir / "8abc_B.develop291.domain_summary.xml").write_text("<domain/>")
            (domains_dir / "7xyz_C.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            proteins = finder._get_proteins_in_batch("ecod_batch_001_test")

            assert len(proteins) == 3
            assert "8ovp_A" in proteins
            assert "8abc_B" in proteins
            assert "7xyz_C" in proteins

            # Should be sorted
            assert proteins == sorted(proteins)

    def test_get_proteins_in_batch_caching(self):
        """Test that _get_proteins_in_batch() caches results"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()

            (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)

            # First call
            proteins1 = finder._get_proteins_in_batch("ecod_batch_001_test")

            # Second call (should use cache)
            proteins2 = finder._get_proteins_in_batch("ecod_batch_001_test")

            assert proteins1 == proteins2
            assert "ecod_batch_001_test" in finder._batch_cache


@pytest.mark.unit
class TestBatchFinderFindBatch:
    """Test find_batch_for_protein() functionality"""

    def test_find_batch_for_protein_stable_batch(self):
        """Test that stable batches are checked first"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create stable batch
            batch_dir = Path(tmpdir, "ecod_batch_036_20250406_1424")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()
            (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)

            # Should find stable batch
            batch = finder.find_batch_for_protein("8ovp_A", verbose=False)

            assert batch == "ecod_batch_036_20250406_1424"

    def test_find_batch_for_protein_not_found(self):
        """Test find_batch_for_protein() returns None when not found"""
        with tempfile.TemporaryDirectory() as tmpdir:
            finder = BatchFinder(tmpdir)

            batch = finder.find_batch_for_protein("nonexistent_protein", verbose=False)

            assert batch is None

    def test_find_batch_for_protein_single_match(self):
        """Test find_batch_for_protein() with single matching batch"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()
            (domains_dir / "8abc_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            batch = finder.find_batch_for_protein("8abc_A", verbose=False)

            assert batch == "ecod_batch_001_test"

    def test_find_batch_for_protein_multiple_matches(self):
        """Test find_batch_for_protein() chooses most recent with multiple matches"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create two batches with same protein
            for batch_name in ["ecod_batch_001_test", "ecod_batch_002_test"]:
                batch_dir = Path(tmpdir, batch_name)
                batch_dir.mkdir()
                domains_dir = batch_dir / "domains"
                domains_dir.mkdir()
                (domains_dir / "8abc_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            batch = finder.find_batch_for_protein("8abc_A", verbose=False)

            # Should choose most recent (002)
            assert batch == "ecod_batch_002_test"

    def test_find_batch_for_protein_verbose_output(self, capsys):
        """Test find_batch_for_protein() verbose output"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()
            (domains_dir / "8abc_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            batch = finder.find_batch_for_protein("8abc_A", verbose=True)

            captured = capsys.readouterr()
            assert "Searching for 8abc_A" in captured.out
            assert "Found in ecod_batch_001_test" in captured.out


@pytest.mark.unit
class TestBatchFinderSuggestSimilar:
    """Test suggest_similar_proteins() functionality"""

    def test_suggest_similar_proteins_empty(self):
        """Test suggest_similar_proteins() with no matches"""
        with tempfile.TemporaryDirectory() as tmpdir:
            finder = BatchFinder(tmpdir)
            suggestions = finder.suggest_similar_proteins("8xyz_A")

            assert suggestions == []

    def test_suggest_similar_proteins_same_pdb(self):
        """Test suggest_similar_proteins() finds proteins from same PDB"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()

            # Create proteins from same PDB
            (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")
            (domains_dir / "8ovp_B.develop291.domain_summary.xml").write_text("<domain/>")
            (domains_dir / "8ovp_C.develop291.domain_summary.xml").write_text("<domain/>")
            (domains_dir / "8abc_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            suggestions = finder.suggest_similar_proteins("8ovp_D")

            # Should suggest other 8ovp chains
            assert "8ovp_A" in suggestions
            assert "8ovp_B" in suggestions
            assert "8ovp_C" in suggestions
            assert "8ovp_D" not in suggestions  # Don't suggest the query itself
            assert "8abc_A" not in suggestions  # Different PDB

    def test_suggest_similar_proteins_max_limit(self):
        """Test suggest_similar_proteins() respects max_suggestions"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()

            # Create many proteins from same PDB
            for chain in "ABCDEFGHIJ":
                (domains_dir / f"8ovp_{chain}.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            suggestions = finder.suggest_similar_proteins("8ovp_Z", max_suggestions=3)

            # Should limit to 3 suggestions
            assert len(suggestions) <= 3


@pytest.mark.unit
class TestBatchFinderAnalyzeBatches:
    """Test analyze_protein_batches() functionality"""

    def test_analyze_protein_batches_not_found(self):
        """Test analyze_protein_batches() with protein not found"""
        with tempfile.TemporaryDirectory() as tmpdir:
            finder = BatchFinder(tmpdir)
            result = finder.analyze_protein_batches("nonexistent_A")

            assert result["multi_batch"] is False
            assert result["batches"] == []

    def test_analyze_protein_batches_single_batch(self):
        """Test analyze_protein_batches() with protein in one batch"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()
            domains_dir = batch_dir / "domains"
            domains_dir.mkdir()
            (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            result = finder.analyze_protein_batches("8ovp_A")

            assert result["multi_batch"] is False
            assert len(result["batches"]) == 1
            assert "ecod_batch_001_test" in result["batches"]

    def test_analyze_protein_batches_multiple_batches(self):
        """Test analyze_protein_batches() with protein in multiple batches"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create protein in multiple batches
            for batch_name in ["ecod_batch_001_test", "ecod_batch_002_test", "ecod_batch_003_test"]:
                batch_dir = Path(tmpdir, batch_name)
                batch_dir.mkdir()
                domains_dir = batch_dir / "domains"
                domains_dir.mkdir()
                (domains_dir / "8ovp_A.develop291.domain_summary.xml").write_text("<domain/>")

            finder = BatchFinder(tmpdir)
            result = finder.analyze_protein_batches("8ovp_A")

            assert result["multi_batch"] is True
            assert len(result["batches"]) == 3
            assert "ecod_batch_001_test" in result["batches"]
            assert "ecod_batch_002_test" in result["batches"]
            assert "ecod_batch_003_test" in result["batches"]


@pytest.mark.unit
class TestPyEcodMiniConfigInitialization:
    """Test PyEcodMiniConfig initialization"""

    def test_config_initialization(self):
        """Test that PyEcodMiniConfig can be initialized"""
        config = PyEcodMiniConfig()

        assert config.base_dir == Path("/data/ecod/pdb_updates/batches")
        assert config.test_data_dir.name == "test_data"
        assert config.output_dir == Path("/tmp")

        # Should have batch finder
        assert isinstance(config.batch_finder, BatchFinder)

    def test_config_reference_files(self):
        """Test that config defines reference file paths"""
        config = PyEcodMiniConfig()

        assert config.domain_lengths_file.name == "domain_lengths.csv"
        assert config.protein_lengths_file.name == "protein_lengths.csv"
        assert config.domain_definitions_file.name == "domain_definitions.csv"
        assert config.reference_blacklist_file.name == "reference_blacklist.csv"


@pytest.mark.unit
class TestPyEcodMiniConfigBatchResolution:
    """Test batch name resolution"""

    def test_resolve_batch_name_full_name(self):
        """Test _resolve_batch_name() with full batch name"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_036_20250406_1424")
            batch_dir.mkdir()

            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            resolved = config._resolve_batch_name("ecod_batch_036_20250406_1424")

            assert resolved == "ecod_batch_036_20250406_1424"

    def test_resolve_batch_name_number(self):
        """Test _resolve_batch_name() with batch number"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_036_20250406_1424")
            batch_dir.mkdir()

            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            resolved = config._resolve_batch_name("036")

            assert resolved == "ecod_batch_036_20250406_1424"

    def test_resolve_batch_name_number_with_padding(self):
        """Test _resolve_batch_name() pads batch numbers correctly"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_20250101_0000")
            batch_dir.mkdir()

            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            # Test with "1" -> should pad to "001"
            resolved = config._resolve_batch_name("1")

            assert resolved == "ecod_batch_001_20250101_0000"

    def test_resolve_batch_name_not_found(self):
        """Test _resolve_batch_name() raises error when not found"""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            with pytest.raises(ValueError) as exc_info:
                config._resolve_batch_name("999")

            assert "No batch found" in str(exc_info.value)


@pytest.mark.unit
class TestPyEcodMiniConfigGetBatch:
    """Test get_batch_for_protein() functionality"""

    def test_get_batch_for_protein_with_explicit_batch(self):
        """Test get_batch_for_protein() with explicit batch_id"""
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir, "ecod_batch_001_test")
            batch_dir.mkdir()

            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            batch = config.get_batch_for_protein("8ovp_A", batch_id="ecod_batch_001_test")

            assert batch == "ecod_batch_001_test"

    def test_get_batch_for_protein_not_found_error(self):
        """Test get_batch_for_protein() raises helpful error when not found"""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PyEcodMiniConfig()
            config.base_dir = Path(tmpdir)

            with pytest.raises(FileNotFoundError) as exc_info:
                config.get_batch_for_protein("nonexistent_A")

            error_msg = str(exc_info.value)
            assert "not found in any batch" in error_msg


@pytest.mark.integration
class TestPyEcodMiniConfigIntegration:
    """Integration tests for PyEcodMiniConfig with real data"""

    def test_config_list_available_batches(self):
        """Test list_available_batches() with real data"""
        base_dir = "/data/ecod/pdb_updates/batches"

        if not Path(base_dir).exists():
            pytest.skip("ECOD batch directory not available")

        config = PyEcodMiniConfig()
        batch_info = config.list_available_batches()

        # Should return list of tuples (batch_name, protein_count)
        assert isinstance(batch_info, list)

        if batch_info:
            assert isinstance(batch_info[0], tuple)
            assert len(batch_info[0]) == 2
            batch_name, protein_count = batch_info[0]
            assert isinstance(batch_name, str)
            assert isinstance(protein_count, int)

    def test_config_validate_setup(self):
        """Test validate_setup() checks configuration"""
        config = PyEcodMiniConfig()

        is_valid, issues = config.validate_setup()

        # Should return bool and list of issues
        assert isinstance(is_valid, bool)
        assert isinstance(issues, list)

        # If invalid, should have issues
        if not is_valid:
            assert len(issues) > 0

    def test_config_get_paths_for_protein(self):
        """Test get_paths_for_protein() returns correct paths"""
        base_dir = "/data/ecod/pdb_updates/batches"

        if not Path(base_dir).exists():
            pytest.skip("ECOD batch directory not available")

        config = PyEcodMiniConfig()

        try:
            paths = config.get_paths_for_protein("8ovp_A", batch_id="036")
        except (FileNotFoundError, ValueError):
            pytest.skip("Batch 036 or protein 8ovp_A not available")

        # Check returned paths
        assert "batch_dir" in paths
        assert "batch_name" in paths
        assert "domain_summary" in paths
        assert "blast_xml" in paths
        assert "output" in paths

        # Paths should be Path objects
        assert isinstance(paths["batch_dir"], Path)
        assert isinstance(paths["domain_summary"], Path)
