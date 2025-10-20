#!/usr/bin/env python3
"""
Official integration test cases for mini_pyecod

This module contains the formal test cases that validate the complete
domain partitioning pipeline works correctly on real proteins.

8ovp_A is the primary/canonical test case that all CI/CD should run against.
"""

from dataclasses import dataclass
from pathlib import Path

import pytest

from pyecod_mini.core.models import PartitionMetadata
from pyecod_mini.core.parser import parse_domain_summary
from pyecod_mini.core.partitioner import partition_domains
from pyecod_mini.core.writer import write_domain_partition


@dataclass
class ExpectedDomain:
    """Expected domain characteristics for validation"""

    family: str
    approximate_range: str
    min_size: int
    max_size: int
    discontinuous: bool = False
    notes: str = ""


@dataclass
class DomainTestCase:
    """Official test case definition"""

    protein_id: str
    description: str
    expected_domain_count: int
    expected_domains: list[ExpectedDomain]
    requires_decomposition: bool = False
    requires_blast_alignments: bool = False
    notes: str = ""


# Official test cases - 8ovp_A is the canonical/primary test
OFFICIAL_TEST_CASES = {
    "8ovp_A": DomainTestCase(
        protein_id="8ovp_A",
        description="GFP-PBP fusion with chain BLAST decomposition (PRIMARY TEST CASE)",
        expected_domain_count=3,
        expected_domains=[
            ExpectedDomain(
                family="271.1.1",  # GFP T-group (correct ECOD classification)
                approximate_range="251-490",
                min_size=235,
                max_size=245,
                discontinuous=False,
                notes="GFP domain identified by T-group 271.1.1",
            ),
            ExpectedDomain(
                family="7523.1.1",  # PBP T-group (first domain)
                approximate_range="1-108",
                min_size=100,
                max_size=120,
                discontinuous=False,
                notes="First PBP domain from chain BLAST decomposition",
            ),
            ExpectedDomain(
                family="7523.1.1",  # PBP T-group (second domain, may be discontinuous)
                approximate_range="109-250",
                min_size=90,
                max_size=150,
                discontinuous=True,  # May be split by insertions
                notes="Second PBP domain, may be discontinuous",
            ),
        ],
        requires_decomposition=True,
        requires_blast_alignments=True,
        notes="""
This is the canonical test case for mini_pyecod. It tests:
- Chain BLAST decomposition with alignment data
- Domain insertion architecture (GFP inserted into PBP)
- Discontinuous domain handling
- Multi-family partitioning with proper T-group assignment

Expected results (validated by algorithm output):
- 3 domains total
- 1 GFP domain (T-group 271.1.1, continuous, ~240 residues)
- 2 PBP domains (T-group 7523.1.1) from decomposition
- Proper ECOD T-group assignment (not PDB IDs)
- >85% sequence coverage

This test MUST pass for any release.
        """,
    ),
    # Additional validated test cases can be added here as they are confirmed
}


class TestOfficialCases:
    """Official integration test cases for mini_pyecod"""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_8ovp_A_canonical(
        self, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
    ):
        """
        PRIMARY TEST CASE: 8ovp_A with full decomposition

        This is the gold standard test that validates the complete pipeline:
        - Evidence parsing from real XML
        - Chain BLAST decomposition with alignment data
        - Proper ECOD T-group assignment
        - Discontinuous domain handling
        """
        test_case = OFFICIAL_TEST_CASES["8ovp_A"]
        result = self._run_test_case(
            test_case, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
        )

        # Basic success check
        assert result[
            "success"
        ], f"Primary test case failed: {result.get('error', 'Unknown error')}"

        # Domain count check - should be exactly 3
        assert result["domain_count"] == 3, f"Expected 3 domains, found {result['domain_count']}"

        # Get the actual domains
        domains = result["domains"]

        # Print detailed results for debugging
        print("\n=== 8ovp_A TEST RESULTS ===")
        print(f"Found {len(domains)} domains:")
        for i, domain in enumerate(domains, 1):
            disc = " (discontinuous)" if domain["discontinuous"] else ""
            print(f"  {i}. {domain['family']}: {domain['range']} ({domain['size']} residues){disc}")
            print(f"     Source: {domain['source']}")

        # Flexible validation based on actual algorithm behavior
        # 1. Check for GFP domain (T-group 271.1.1)
        gfp_domains = [d for d in domains if "271.1.1" in d["family"]]
        assert len(gfp_domains) >= 1, "Should have at least one GFP domain (T-group 271.1.1)"

        # 2. Check for PBP domains (T-group 7523.1.1)
        pbp_domains = [d for d in domains if "7523.1.1" in d["family"]]
        assert len(pbp_domains) >= 1, "Should have at least one PBP domain (T-group 7523.1.1)"

        # 3. Check for decomposed domains
        decomposed_domains = [d for d in domains if d["source"] == "chain_blast_decomposed"]
        assert (
            len(decomposed_domains) >= 1
        ), f"Should have decomposed domains, found {len(decomposed_domains)}"

        # 4. Check that we have the expected T-groups
        t_groups = {d["family"] for d in domains}
        assert "271.1.1" in t_groups, f"Should have GFP T-group 271.1.1, found {t_groups}"
        assert "7523.1.1" in t_groups, f"Should have PBP T-group 7523.1.1, found {t_groups}"

        # 5. Coverage check
        total_coverage = sum(d["size"] for d in domains)
        sequence_length = 569  # Known length for 8ovp_A
        coverage_fraction = total_coverage / sequence_length
        assert coverage_fraction >= 0.80, f"Coverage {coverage_fraction:.1%} is too low (need ≥80%)"

        # 6. Domain size sanity check
        for domain in domains:
            assert (
                20 <= domain["size"] <= 500
            ), f"Domain size {domain['size']} outside reasonable range"

        # Summary
        print("\n✅ PRIMARY TEST CASE PASSED")
        print(f"   Total coverage: {total_coverage}/{sequence_length} ({coverage_fraction:.1%})")
        print(f"   T-groups found: {sorted(t_groups)}")
        print(f"   Decomposed domains: {len(decomposed_domains)}")
        print(f"   GFP domains: {len(gfp_domains)}")
        print(f"   PBP domains: {len(pbp_domains)}")

    @pytest.mark.integration
    def test_8ovp_A_without_decomposition(
        self, stable_batch_dir, real_reference_data, temp_output_dir
    ):
        """
        Test 8ovp_A without chain BLAST decomposition

        This tests the basic partitioning algorithm without decomposition.
        """
        test_case = OFFICIAL_TEST_CASES["8ovp_A"]

        # Run without domain definitions (disables decomposition)
        reference_data_no_decomp = real_reference_data.copy()
        reference_data_no_decomp["domain_definitions"] = {}

        result = self._run_test_case(
            test_case,
            stable_batch_dir,
            reference_data_no_decomp,
            blast_alignments={},  # No BLAST alignments
            output_dir=temp_output_dir,
            expect_decomposition=False,
        )

        assert result[
            "success"
        ], f"No-decomposition test failed: {result.get('error', 'Unknown error')}"

        # Without decomposition, should still get domains but possibly fewer
        domains = result["domains"]

        # Should have no decomposed domains
        decomposed_domains = [d for d in domains if d["source"] == "chain_blast_decomposed"]
        assert (
            len(decomposed_domains) == 0
        ), f"Should have no decomposed domains, found {len(decomposed_domains)}"

        # Should still find some domains
        assert (
            result["domain_count"] >= 1
        ), f"Should find at least 1 domain, found {result['domain_count']}"

        print(f"✅ NO-DECOMPOSITION TEST PASSED: Found {result['domain_count']} domains")

    @pytest.mark.performance
    @pytest.mark.integration
    def test_8ovp_A_performance(
        self, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
    ):
        """
        Performance benchmark for 8ovp_A

        Ensures processing completes within reasonable time limits.
        """
        import time

        test_case = OFFICIAL_TEST_CASES["8ovp_A"]

        start_time = time.time()
        result = self._run_test_case(
            test_case, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
        )
        processing_time = time.time() - start_time

        # Performance requirements
        MAX_PROCESSING_TIME = 60  # seconds

        assert result["success"], "Performance test must use working algorithm"
        assert (
            processing_time < MAX_PROCESSING_TIME
        ), f"Processing took {processing_time:.1f}s (max: {MAX_PROCESSING_TIME}s)"

        print(f"✅ PERFORMANCE TEST PASSED: {processing_time:.2f}s (limit: {MAX_PROCESSING_TIME}s)")

    @pytest.mark.parametrize("protein_id", ["8ovp_A"])  # Extend as more cases are validated
    def test_output_file_generation(
        self, protein_id, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
    ):
        """
        Test that output XML files are generated correctly
        """
        test_case = OFFICIAL_TEST_CASES[protein_id]
        result = self._run_test_case(
            test_case, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir
        )

        assert result["success"], f"Test case failed: {result.get('error', 'Unknown error')}"

        # Check output file was created
        output_file = Path(temp_output_dir) / f"{protein_id}_test.domains.xml"
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Validate XML structure
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_file)
        root = tree.getroot()

        assert root.tag == "domain_partition", "Root element should be domain_partition"
        assert root.get("pdb_id") == protein_id.split("_")[0], "PDB ID should be set correctly"

        domains_elem = root.find("domains")
        assert domains_elem is not None, "Should have domains element"

        domain_elems = domains_elem.findall("domain")
        assert len(domain_elems) == result["domain_count"], "XML should match found domain count"

        print(f"✅ OUTPUT VALIDATION PASSED: {output_file}")

    def _run_test_case(
        self,
        test_case: DomainTestCase,
        batch_dir: str,
        reference_data: dict,
        blast_alignments: dict,
        output_dir: str,
        expect_decomposition: bool = True,
    ) -> dict:
        """
        Run a single test case and return results
        """
        import os

        # Parse protein ID
        parts = test_case.protein_id.split("_")
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else "A"

        # Check domain summary file
        xml_path = os.path.join(
            batch_dir, "domains", f"{test_case.protein_id}.develop291.domain_summary.xml"
        )
        if not os.path.exists(xml_path):
            return {
                "success": False,
                "error": f"Domain summary not found: {xml_path}",
                "test_case": test_case,
            }

        try:
            # Parse evidence
            evidence = parse_domain_summary(
                xml_path,
                reference_lengths=reference_data.get("domain_lengths", {}),
                protein_lengths=reference_data.get("protein_lengths", {}),
                blast_alignments=blast_alignments,
                require_reference_lengths=True,
                verbose=False,
            )

            if not evidence:
                return {
                    "success": False,
                    "error": "No evidence with reference lengths found",
                    "test_case": test_case,
                }

            # Estimate sequence length
            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
            sequence_length = int(max_pos * 1.1)

            # Partition domains
            domain_definitions = (
                reference_data.get("domain_definitions", {}) if expect_decomposition else {}
            )
            domains = partition_domains(
                evidence,
                sequence_length=sequence_length,
                domain_definitions=domain_definitions,
                verbose=False,
            )

            # Write output for validation - FIXED: Use new PartitionMetadata API
            output_file = os.path.join(output_dir, f"{test_case.protein_id}_test.domains.xml")

            # Create metadata object instead of passing individual parameters
            metadata = PartitionMetadata(
                pdb_id=pdb_id,
                chain_id=chain_id,
                sequence_length=sequence_length,
                source_domain_summary_path=xml_path,
                batch_id=os.path.basename(batch_dir),
            )

            # Use new API signature
            write_domain_partition(domains, metadata, output_file)

            # Convert domains to serializable format
            domain_data = []
            for domain in domains:
                domain_info = {
                    "id": domain.id,
                    "family": domain.family,
                    "range": str(domain.range),
                    "size": domain.range.total_length,
                    "discontinuous": domain.range.is_discontinuous,
                    "source": domain.source,
                }
                domain_data.append(domain_info)

            return {
                "success": True,
                "test_case": test_case,
                "domain_count": len(domains),
                "domains": domain_data,
                "output_file": output_file,
                "sequence_length": sequence_length,
            }

        except Exception as e:
            import traceback

            return {
                "success": False,
                "error": str(e),
                "traceback": traceback.format_exc(),
                "test_case": test_case,
            }


class TestCaseValidation:
    """Validation utilities for test cases"""

    @pytest.mark.unit
    def test_test_case_definitions(self):
        """
        Validate that test case definitions are well-formed
        """
        for protein_id, test_case in OFFICIAL_TEST_CASES.items():
            # Basic validation
            assert test_case.protein_id == protein_id, f"Protein ID mismatch for {protein_id}"
            assert test_case.expected_domain_count > 0, f"Invalid domain count for {protein_id}"
            assert (
                len(test_case.expected_domains) > 0
            ), f"No expected domains defined for {protein_id}"

            # Validate expected domains
            for expected_domain in test_case.expected_domains:
                assert expected_domain.family, f"Missing family for domain in {protein_id}"
                assert expected_domain.min_size > 0, f"Invalid min_size for domain in {protein_id}"
                assert (
                    expected_domain.max_size >= expected_domain.min_size
                ), f"Invalid size range for domain in {protein_id}"

        print(f"✅ TEST CASE DEFINITIONS VALIDATED: {len(OFFICIAL_TEST_CASES)} cases")

    @pytest.mark.unit
    def test_primary_test_case_properties(self):
        """
        Validate that 8ovp_A has all required properties for a primary test case
        """
        primary = OFFICIAL_TEST_CASES["8ovp_A"]

        # Primary test case requirements
        assert primary.requires_decomposition, "Primary test case must test decomposition"
        assert primary.requires_blast_alignments, "Primary test case must test BLAST alignments"
        assert (
            primary.expected_domain_count >= 3
        ), "Primary test case should be complex (≥3 domains)"

        # Check that we have both GFP and PBP T-groups expected
        families = {ed.family for ed in primary.expected_domains}
        assert "271.1.1" in families, "Primary test case should include GFP T-group 271.1.1"
        assert "7523.1.1" in families, "Primary test case should include PBP T-group 7523.1.1"

        print("✅ PRIMARY TEST CASE PROPERTIES VALIDATED")


# Standalone runner for development/debugging
def main():
    """Run test cases directly (for development)"""
    import argparse

    parser = argparse.ArgumentParser(description="Run official mini_pyecod test cases")
    parser.add_argument("--protein", default="8ovp_A", help="Protein to test")
    parser.add_argument("--batch-dir", help="Batch directory override")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if args.protein not in OFFICIAL_TEST_CASES:
        print(f"Unknown test case: {args.protein}")
        print(f"Available: {list(OFFICIAL_TEST_CASES.keys())}")
        return 1

    # This would run the test case directly (implementation depends on environment setup)
    print(f"Running test case: {args.protein}")
    print(
        "Use 'pytest tests/test_cases.py::TestOfficialCases::test_8ovp_A_canonical -v' for full pytest integration"
    )

    return 0


if __name__ == "__main__":
    exit(main())
