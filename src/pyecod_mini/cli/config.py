#!/usr/bin/env python3
"""
Configuration and batch finding for pyecod_mini CLI

Handles smart batch detection, path resolution, and configuration management.
"""

from pathlib import Path
from typing import Optional


class BatchFinder:
    """Smart batch finder for proteins"""

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self._batch_cache = {}

        # Known stable batches for test cases
        self.stable_batches = {
            "8ovp_A": "ecod_batch_036_20250406_1424",  # Our validated test case
        }

    def find_batch_for_protein(self, protein_id: str, verbose: bool = False) -> Optional[str]:
        """Find which batch contains the protein"""

        # Check for known stable test cases first
        if protein_id in self.stable_batches:
            stable_batch = self.stable_batches[protein_id]
            if self._protein_exists_in_batch(protein_id, stable_batch):
                if verbose:
                    print(f"Using stable batch for {protein_id}: {stable_batch}")
                return stable_batch
            if verbose:
                print(f"WARNING: Stable batch {stable_batch} doesn't contain {protein_id}")

        # Search all available batches
        available_batches = self._get_available_batches()
        found_batches = []

        if verbose:
            print(f"Searching for {protein_id} across {len(available_batches)} batches...")

        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)
                if verbose:
                    print(f"  ✓ Found in {batch_name}")

        if not found_batches:
            if verbose:
                print(f"  ✗ {protein_id} not found in any batch")
            return None

        if len(found_batches) == 1:
            return found_batches[0]

        # Multiple batches found - warn and use most recent
        if verbose:
            print(f"  ⚠️  Found in {len(found_batches)} batches: {found_batches}")
            print(f"  → Using most recent: {found_batches[-1]}")
            print("  💡 Use --batch-id to specify a particular batch")

        return found_batches[-1]  # Most recent

    def suggest_similar_proteins(self, protein_id: str, max_suggestions: int = 5) -> list[str]:
        """Suggest similar protein IDs that exist"""
        pdb_id = protein_id.split("_")[0] if "_" in protein_id else protein_id[:4]

        suggestions = []
        for batch_name in self._get_available_batches():
            proteins = self._get_proteins_in_batch(batch_name)
            same_pdb = [p for p in proteins if p.startswith(pdb_id)]
            suggestions.extend(same_pdb)

            if len(suggestions) >= max_suggestions:
                break

        unique_suggestions = list(dict.fromkeys(suggestions))[:max_suggestions]
        return [s for s in unique_suggestions if s != protein_id]

    def analyze_protein_batches(self, protein_id: str) -> dict[str, any]:
        """Analyze a protein across multiple batches"""
        available_batches = self._get_available_batches()
        found_batches = []

        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)

        if len(found_batches) <= 1:
            return {"multi_batch": False, "batches": found_batches}

        return {"multi_batch": True, "batches": found_batches}

    def _get_available_batches(self) -> list[str]:
        """Get list of available batch directories"""
        if not self.base_dir.exists():
            return []

        # Include both ecod_batch_* AND alt_rep_batch_* directories
        batch_dirs = [
            d.name
            for d in self.base_dir.iterdir()
            if d.is_dir()
            and (d.name.startswith("ecod_batch_") or d.name.startswith("alt_rep_batch_"))
        ]

        return sorted(batch_dirs)

    def _protein_exists_in_batch(self, protein_id: str, batch_name: str) -> bool:
        """Check if protein exists in a specific batch"""
        batch_dir = self.base_dir / batch_name
        domain_file = batch_dir / "domains" / f"{protein_id}.develop291.domain_summary.xml"
        return domain_file.exists()

    def _get_proteins_in_batch(self, batch_name: str) -> list[str]:
        """Get list of proteins in a batch (cached)"""
        if batch_name in self._batch_cache:
            return self._batch_cache[batch_name]

        batch_dir = self.base_dir / batch_name
        domains_dir = batch_dir / "domains"

        if not domains_dir.exists():
            self._batch_cache[batch_name] = []
            return []

        proteins = []
        for domain_file in domains_dir.glob("*.develop291.domain_summary.xml"):
            protein_id = domain_file.stem.replace(".develop291.domain_summary", "")
            proteins.append(protein_id)

        self._batch_cache[batch_name] = sorted(proteins)
        return self._batch_cache[batch_name]


class PyEcodMiniConfig:
    """Configuration manager with integrated batch detection"""

    def __init__(self):
        # Try to find project root (look for pyproject.toml)
        current_dir = Path(__file__).parent
        project_root = None
        for parent in [current_dir] + list(current_dir.parents):
            if (parent / "pyproject.toml").exists():
                project_root = parent
                break

        if project_root is None:
            project_root = Path(__file__).parent.parent.parent

        self.base_dir = Path("/data/ecod/pdb_updates/batches")
        self.test_data_dir = project_root / "test_data"
        self.output_dir = Path("/tmp")

        # Default reference files
        self.domain_lengths_file = self.test_data_dir / "domain_lengths.csv"
        self.protein_lengths_file = self.test_data_dir / "protein_lengths.csv"
        self.domain_definitions_file = self.test_data_dir / "domain_definitions.csv"
        self.reference_blacklist_file = self.test_data_dir / "reference_blacklist.csv"

        # Batch finder
        self.batch_finder = BatchFinder(str(self.base_dir))

        # Visualization settings
        self.pdb_repo_path = "/usr2/pdb/data"
        self.visualization_output_dir = "/tmp/pymol_comparison"

    def get_batch_for_protein(
        self, protein_id: str, batch_id: Optional[str] = None, verbose: bool = False
    ) -> str:
        """Get the right batch for a protein"""
        if batch_id is not None:
            # Explicit batch specified - validate and resolve
            return self._resolve_batch_name(batch_id)

        # Smart detection
        found_batch = self.batch_finder.find_batch_for_protein(protein_id, verbose)

        if found_batch is None:
            # Provide helpful error message
            suggestions = self.batch_finder.suggest_similar_proteins(protein_id)
            error_msg = f"Protein {protein_id} not found in any batch"

            if suggestions:
                error_msg += f"\nSimilar proteins available: {suggestions[:3]}"

            available_batches = self.batch_finder._get_available_batches()
            if available_batches:
                error_msg += f"\nAvailable batches: {available_batches[-3:]}"
                error_msg += f"\nTo specify a batch: pyecod-mini {protein_id} --batch-id BATCH_NAME"

            raise FileNotFoundError(error_msg)

        return found_batch

    def _resolve_batch_name(self, batch_id: str) -> str:
        """Convert batch_id to full batch name"""
        if batch_id.isdigit():
            # Number like "036" -> find "ecod_batch_036_*"
            pattern = f"ecod_batch_{batch_id.zfill(3)}_*"
            matches = list(self.base_dir.glob(pattern))
            if matches:
                return matches[0].name
            msg = f"No batch found matching pattern: {pattern}"
            raise ValueError(msg)
        # Assume it's already a full batch name
        if (self.base_dir / batch_id).exists():
            return batch_id
        msg = f"Batch directory not found: {batch_id}"
        raise ValueError(msg)

    def get_batch_dir(self, batch_name: str) -> Path:
        """Get batch directory path from batch name"""
        return self.base_dir / batch_name

    def get_paths_for_protein(
        self, protein_id: str, batch_id: Optional[str] = None, verbose: bool = False
    ) -> dict[str, Path]:
        """Get all file paths for a protein"""

        # Find the right batch for this protein
        batch_name = self.get_batch_for_protein(protein_id, batch_id, verbose)
        batch_dir = self.get_batch_dir(batch_name)

        if verbose:
            print(f"Using batch: {batch_name}")

        return {
            "batch_dir": batch_dir,
            "batch_name": batch_name,
            "domain_summary": batch_dir / "domains" / f"{protein_id}.develop291.domain_summary.xml",
            "blast_xml": batch_dir / "blast" / "chain" / f"{protein_id}.develop291.xml",
            "blast_dir": batch_dir / "blast" / "chain",
            "domain_lengths": self.domain_lengths_file,
            "protein_lengths": self.protein_lengths_file,
            "domain_definitions": self.domain_definitions_file,
            "output": self.output_dir / f"{protein_id}_mini.domains.xml",
            "old_domains": batch_dir / "domains" / f"{protein_id}.develop291.domains.xml",
        }

    def list_available_batches(self) -> list[tuple[str, int]]:
        """List all available batch directories with protein counts"""
        batches = self.batch_finder._get_available_batches()

        batch_info = []
        for batch_name in batches:
            proteins = self.batch_finder._get_proteins_in_batch(batch_name)
            batch_info.append((batch_name, len(proteins)))

        return batch_info

    def validate_setup(self, verbose: bool = False) -> tuple[bool, list[str]]:
        """Validate that the configuration is usable"""
        issues = []

        # Check base directory
        if not self.base_dir.exists():
            issues.append(f"Base directory not found: {self.base_dir}")
            return False, issues

        # Check if any batches exist
        available_batches = self.batch_finder._get_available_batches()
        if not available_batches:
            issues.append("No batch directories found")

        # Check test data files
        for name, path in [
            ("domain lengths", self.domain_lengths_file),
            ("protein lengths", self.protein_lengths_file),
        ]:
            if not path.exists():
                issues.append(f"{name.title()} file not found: {path}")

        # Domain definitions are optional but recommended
        if not self.domain_definitions_file.exists():
            issues.append(
                f"Domain definitions file not found (chain BLAST decomposition disabled): {self.domain_definitions_file}"
            )

        # Blacklist is optional
        if not self.reference_blacklist_file.exists() and verbose:
            print(f"No reference blacklist found: {self.reference_blacklist_file}")

        return len(issues) == 0, issues
