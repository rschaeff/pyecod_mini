# Production Framework Design

**Date**: 2025-10-19
**System**: pyECOD Mini Production Processing
**Version**: 2.0.0

## Overview

This document details the design of the production framework for processing ~40,000 representative proteins using the validated pyECOD Mini algorithm with SLURM cluster integration.

## Design Principles

Following lessons from CLAUDE.md:

1. ✅ **Simple Over Complex**: No service-oriented architecture, no complex orchestration
2. ✅ **Algorithm First**: Build around proven mini algorithm, not before
3. ✅ **Proven Technologies**: SLURM (battle-tested), PostgreSQL (reliable), filesystem-based tracking
4. ✅ **Observable**: Real-time monitoring, comprehensive logging, quality metrics
5. ✅ **Maintainable**: Clean code, clear separation of concerns, good documentation

## System Architecture

### High-Level Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                        Production Controller                         │
│  (Python CLI: pyecod-mini production ...)                           │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                    ┌──────────────┼──────────────┐
                    │              │              │
                    v              v              v
         ┌─────────────┐  ┌──────────────┐  ┌──────────┐
         │   Scanner   │  │ SLURM Manager│  │ Monitor  │
         └─────────────┘  └──────────────┘  └──────────┘
                 │                │                │
                 v                v                v
         ┌──────────────────────────────────────────────┐
         │      SQLite Tracking Database                │
         │  (Job status, metadata, quality metrics)     │
         └──────────────────────────────────────────────┘
                                   │
                                   v
         ┌──────────────────────────────────────────────┐
         │        SLURM Cluster (50 concurrent jobs)    │
         │                                               │
         │  ┌─────────┐ ┌─────────┐      ┌─────────┐  │
         │  │ Worker  │ │ Worker  │ ...  │ Worker  │  │
         │  │  Job 1  │ │  Job 2  │      │  Job N  │  │
         │  └─────────┘ └─────────┘      └─────────┘  │
         │       │            │                │        │
         └───────┼────────────┼────────────────┼────────┘
                 │            │                │
                 v            v                v
         ┌──────────────────────────────────────────────┐
         │     pyecod-mini core algorithm               │
         │  (Domain partitioning, boundary opt)         │
         └──────────────────────────────────────────────┘
                                   │
                    ┌──────────────┴──────────────┐
                    │                             │
                    v                             v
         ┌───────────────────┐         ┌──────────────────┐
         │  XML Output Files │         │  Tracking Update │
         │  (batch/mini/)    │         │  (SQLite DB)     │
         └───────────────────┘         └──────────────────┘
                    │
                    v
         ┌──────────────────────────────────────────────┐
         │           Database Importer                   │
         │  (PostgreSQL: pdb_analysis schema)           │
         └──────────────────────────────────────────────┘
                                   │
                                   v
         ┌──────────────────────────────────────────────┐
         │    PostgreSQL Production Database            │
         │  - partition_proteins                        │
         │  - partition_domains                         │
         │  - partition_evidence                        │
         └──────────────────────────────────────────────┘
```

## Component Specifications

### 1. Scanner Component

**Purpose**: Identify proteins to process, filter by criteria, check existing results

**Module**: `src/pyecod_mini/production/scanner.py`

**Key Classes**:

```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional, Set
from datetime import datetime

@dataclass
class ProteinTask:
    """Represents a single protein processing task"""
    protein_id: str
    pdb_id: str
    chain_id: str
    batch_name: str
    batch_dir: Path

    # Input files
    domain_summary_path: Path
    blast_xml_path: Optional[Path] = None

    # Status tracking
    is_representative: bool = False
    existing_output: Optional[Path] = None
    previous_attempt: Optional[datetime] = None
    retry_count: int = 0

    # Priority (higher = process first)
    priority: int = 0

@dataclass
class ScanResult:
    """Result of scanning batch directories"""
    total_proteins_found: int
    representative_count: int
    tasks_to_process: List[ProteinTask]
    tasks_skipped: List[ProteinTask]  # Already processed
    tasks_failed: List[ProteinTask]   # Failed previously

    scan_timestamp: datetime = field(default_factory=datetime.now)
    scan_duration_seconds: float = 0.0

    def summary(self) -> str:
        """Human-readable summary"""
        return f"""
Scan Results:
  Total proteins: {self.total_proteins_found}
  Representatives: {self.representative_count}
  To process: {len(self.tasks_to_process)}
  Skipped (done): {len(self.tasks_skipped)}
  Failed (retry): {len(self.tasks_failed)}
        """

class BatchScanner:
    """Scan batch directories and identify proteins to process"""

    def __init__(
        self,
        batch_base_dir: Path,
        db_connection: Optional[DatabaseConnection] = None,
        tracking_db: Optional[TrackingDatabase] = None
    ):
        self.batch_base_dir = batch_base_dir
        self.db_connection = db_connection
        self.tracking_db = tracking_db

    def scan_batches(
        self,
        batch_patterns: List[str] = ["ecod_batch_*"],
        reps_only: bool = True,
        skip_existing: bool = True,
        max_proteins: Optional[int] = None
    ) -> ScanResult:
        """
        Scan batch directories and identify proteins to process

        Args:
            batch_patterns: Glob patterns for batch directories
            reps_only: Only include representative proteins
            skip_existing: Skip proteins with existing output
            max_proteins: Limit number of proteins (for testing)

        Returns:
            ScanResult with identified tasks
        """
        ...

    def get_representative_proteins(self) -> Set[str]:
        """Query database for representative protein IDs"""
        if not self.db_connection:
            return set()

        # SQL: SELECT source_id FROM ecod_schema.process_status WHERE is_representative = TRUE
        ...

    def check_existing_output(self, protein_id: str, batch_name: str) -> Optional[Path]:
        """Check if output already exists for protein"""
        output_path = self.batch_base_dir / batch_name / "mini_domains" / f"{protein_id}.mini.domains.xml"
        return output_path if output_path.exists() else None

    def prioritize_tasks(self, tasks: List[ProteinTask]) -> List[ProteinTask]:
        """
        Prioritize tasks by various criteria

        Priority factors:
        1. Failed tasks (for retry)
        2. High-quality evidence (from tracking DB)
        3. Simple cases (fewer domains expected)
        4. Batch order
        """
        ...
```

**Usage Example**:
```python
scanner = BatchScanner(
    batch_base_dir=Path("/data/ecod/pdb_updates/batches"),
    db_connection=db,
    tracking_db=tracking
)

result = scanner.scan_batches(
    batch_patterns=["ecod_batch_0*"],
    reps_only=True,
    skip_existing=True,
    max_proteins=1000  # For testing
)

print(result.summary())
# Scan Results:
#   Total proteins: 42,384
#   Representatives: 39,127
#   To process: 1,000
#   Skipped (done): 38,127
#   Failed (retry): 0
```

### 2. SLURM Manager

**Purpose**: Submit, monitor, and manage SLURM jobs for parallel processing

**Module**: `src/pyecod_mini/production/slurm_manager.py`

**Key Classes**:

```python
from enum import Enum

class JobStatus(Enum):
    """SLURM job status"""
    PENDING = "pending"      # Submitted, waiting in queue
    RUNNING = "running"      # Currently executing
    COMPLETED = "completed"  # Finished successfully
    FAILED = "failed"        # Exit code != 0
    TIMEOUT = "timeout"      # Exceeded time limit
    CANCELLED = "cancelled"  # User cancelled
    UNKNOWN = "unknown"      # Cannot determine status

@dataclass
class SlurmJob:
    """SLURM job specification and tracking"""
    # Job identification
    job_id: Optional[str] = None  # SLURM job ID (e.g., "12345678")
    protein_id: str = ""
    batch_name: str = ""

    # Job specification
    partition: str = "All"
    time_limit: str = "1:00:00"
    memory: str = "4G"
    cpus: int = 1

    # Paths
    script_path: Optional[Path] = None
    output_log: Optional[Path] = None
    error_log: Optional[Path] = None

    # Status tracking
    status: JobStatus = JobStatus.PENDING
    submit_time: Optional[datetime] = None
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    exit_code: Optional[int] = None

    # Retry tracking
    retry_count: int = 0
    max_retries: int = 2

    def elapsed_time(self) -> Optional[float]:
        """Get elapsed time in seconds"""
        if not self.start_time:
            return None
        end = self.end_time or datetime.now()
        return (end - self.start_time).total_seconds()

    def is_terminal(self) -> bool:
        """Check if job is in terminal state"""
        return self.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]

    def should_retry(self) -> bool:
        """Check if job should be retried"""
        return self.status == JobStatus.FAILED and self.retry_count < self.max_retries

@dataclass
class JobBatch:
    """Collection of related SLURM jobs"""
    name: str
    jobs: List[SlurmJob] = field(default_factory=list)
    created: datetime = field(default_factory=datetime.now)

    def get_stats(self) -> Dict[str, int]:
        """Get statistics for this batch"""
        stats = {}
        for status in JobStatus:
            stats[status.value] = len([j for j in self.jobs if j.status == status])
        return stats

class SlurmManager:
    """Manage SLURM job submission and monitoring"""

    def __init__(
        self,
        config: SlurmConfig,
        tracking_db: TrackingDatabase,
        job_script_dir: Path = Path("/tmp/pyecod_mini_jobs"),
        log_dir: Path = Path("/tmp/pyecod_mini_logs")
    ):
        self.config = config
        self.tracking_db = tracking_db
        self.job_script_dir = job_script_dir
        self.log_dir = log_dir

        # Create directories
        self.job_script_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def submit_job(self, task: ProteinTask) -> SlurmJob:
        """
        Submit a single protein processing job

        Steps:
        1. Generate SLURM script
        2. Submit with sbatch
        3. Record job ID
        4. Update tracking database
        """
        # Generate script
        script_path = self._generate_job_script(task)

        # Submit to SLURM
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise SlurmSubmissionError(f"sbatch failed: {result.stderr}")

        # Parse job ID from output: "Submitted batch job 12345678"
        job_id = result.stdout.strip().split()[-1]

        # Create job record
        job = SlurmJob(
            job_id=job_id,
            protein_id=task.protein_id,
            batch_name=task.batch_name,
            script_path=script_path,
            output_log=self.log_dir / f"{task.protein_id}_{job_id}.out",
            error_log=self.log_dir / f"{task.protein_id}_{job_id}.err",
            submit_time=datetime.now(),
            status=JobStatus.PENDING
        )

        # Record in tracking database
        self.tracking_db.record_job_submission(job)

        return job

    def submit_batch(
        self,
        tasks: List[ProteinTask],
        max_concurrent: int = 50,
        throttle_delay: float = 0.1
    ) -> JobBatch:
        """
        Submit a batch of jobs with rate limiting

        Args:
            tasks: Protein tasks to submit
            max_concurrent: Maximum concurrent jobs
            throttle_delay: Delay between submissions (seconds)
        """
        batch = JobBatch(name=f"batch_{datetime.now():%Y%m%d_%H%M%S}")

        for task in tasks:
            # Check current job count
            while self._count_active_jobs() >= max_concurrent:
                time.sleep(5)  # Wait for jobs to complete

            # Submit job
            job = self.submit_job(task)
            batch.jobs.append(job)

            # Throttle submissions
            time.sleep(throttle_delay)

        return batch

    def monitor_jobs(self, jobs: List[SlurmJob]) -> Dict[str, int]:
        """
        Monitor job status and update records

        Uses: sacct, squeue

        Returns: Status counts
        """
        # Query SLURM for current status
        job_ids = [j.job_id for j in jobs if j.job_id]

        if not job_ids:
            return {}

        # Use sacct to get job status
        result = subprocess.run(
            ["sacct", "-j", ",".join(job_ids),
             "--format=JobID,State,ExitCode,Start,End",
             "--parsable2", "--noheader"],
            capture_output=True,
            text=True
        )

        # Parse output and update job records
        ...

        # Update tracking database
        for job in jobs:
            self.tracking_db.update_job_status(job)

        return self._count_by_status(jobs)

    def retry_failed_jobs(
        self,
        jobs: List[SlurmJob],
        delay: int = 300
    ) -> List[SlurmJob]:
        """
        Retry failed jobs with exponential backoff

        Args:
            jobs: Jobs to check for retry
            delay: Base delay in seconds

        Returns: Newly submitted retry jobs
        """
        retry_jobs = []

        for job in jobs:
            if not job.should_retry():
                continue

            # Exponential backoff
            wait_time = delay * (2 ** job.retry_count)
            time.sleep(wait_time)

            # Resubmit
            task = self._job_to_task(job)
            new_job = self.submit_job(task)
            new_job.retry_count = job.retry_count + 1
            retry_jobs.append(new_job)

        return retry_jobs

    def _generate_job_script(self, task: ProteinTask) -> Path:
        """Generate SLURM job script for a protein"""
        script_path = self.job_script_dir / f"{task.protein_id}.sh"

        script_content = f"""#!/bin/bash
#SBATCH --job-name=mini_{task.protein_id}
#SBATCH --partition={self.config.partition}
#SBATCH --time={self.config.time_limit}
#SBATCH --mem={self.config.memory}
#SBATCH --cpus-per-task={self.config.cpus}
#SBATCH --output={self.log_dir}/%x_%j.out
#SBATCH --error={self.log_dir}/%x_%j.err

# Load environment
source {self.config.env_setup_script}

# Run pyecod-mini
cd {self.config.work_dir}
pyecod-mini {task.protein_id} \\
    --batch-id {task.batch_name} \\
    --output {task.batch_dir}/mini_domains/{task.protein_id}.mini.domains.xml

# Record completion status
echo "Exit code: $?" > {self.log_dir}/{task.protein_id}_$SLURM_JOB_ID.status
"""

        script_path.write_text(script_content)
        script_path.chmod(0o755)

        return script_path

    def _count_active_jobs(self) -> int:
        """Count currently active (pending + running) jobs"""
        result = subprocess.run(
            ["squeue", "-u", os.environ["USER"], "--format=%A", "--noheader"],
            capture_output=True,
            text=True
        )
        return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
```

**Usage Example**:
```python
manager = SlurmManager(
    config=slurm_config,
    tracking_db=tracking,
    job_script_dir=Path("/tmp/pyecod_mini_jobs"),
    log_dir=Path("/tmp/pyecod_mini_logs")
)

# Submit batch of jobs
batch = manager.submit_batch(
    tasks=scan_result.tasks_to_process,
    max_concurrent=50,
    throttle_delay=0.1
)

# Monitor progress
while not all(j.is_terminal() for j in batch.jobs):
    stats = manager.monitor_jobs(batch.jobs)
    print(f"Status: {stats}")
    time.sleep(60)  # Check every minute

# Retry failed jobs
retry_jobs = manager.retry_failed_jobs(batch.jobs)
```

### 3. Tracking Database

**Purpose**: Lightweight SQLite database for job tracking, metadata, quality metrics

**Module**: `src/pyecod_mini/production/tracking_db.py`

**Schema**:

```sql
CREATE TABLE jobs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id TEXT,  -- SLURM job ID
    protein_id TEXT NOT NULL,
    pdb_id TEXT NOT NULL,
    chain_id TEXT NOT NULL,
    batch_name TEXT NOT NULL,

    -- Status
    status TEXT NOT NULL,  -- pending, running, completed, failed
    submit_time TIMESTAMP,
    start_time TIMESTAMP,
    end_time TIMESTAMP,
    exit_code INTEGER,

    -- Retry tracking
    retry_count INTEGER DEFAULT 0,
    parent_job_id TEXT,  -- Original job if this is a retry

    -- Paths
    script_path TEXT,
    output_log TEXT,
    error_log TEXT,
    output_xml TEXT,

    UNIQUE(protein_id, batch_name)
);

CREATE TABLE results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id INTEGER REFERENCES jobs(id),
    protein_id TEXT NOT NULL,

    -- Results
    domain_count INTEGER,
    coverage_fraction REAL,
    total_residues INTEGER,
    assigned_residues INTEGER,

    -- Quality metrics
    avg_confidence REAL,
    evidence_count INTEGER,
    discontinuous_domains INTEGER,

    -- Processing metadata
    processing_time_seconds REAL,
    peak_memory_mb REAL,

    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE quality_flags (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id INTEGER REFERENCES jobs(id),
    protein_id TEXT NOT NULL,

    flag_type TEXT NOT NULL,  -- low_coverage, poor_evidence, etc.
    severity TEXT NOT NULL,   -- warning, error
    message TEXT,

    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_jobs_protein ON jobs(protein_id);
CREATE INDEX idx_jobs_status ON jobs(status);
CREATE INDEX idx_results_protein ON results(protein_id);
```

**Key Class**:

```python
class TrackingDatabase:
    """SQLite database for tracking job status and results"""

    def __init__(self, db_path: Path):
        self.db_path = db_path
        self.conn = sqlite3.connect(str(db_path))
        self._create_schema()

    def record_job_submission(self, job: SlurmJob) -> None:
        """Record job submission"""
        ...

    def update_job_status(self, job: SlurmJob) -> None:
        """Update job status"""
        ...

    def record_result(self, protein_id: str, result: PartitionResult) -> None:
        """Record processing result and quality metrics"""
        ...

    def get_stats(self) -> Dict[str, Any]:
        """Get overall processing statistics"""
        ...

    def get_failed_jobs(self) -> List[SlurmJob]:
        """Get all failed jobs for retry"""
        ...
```

### 4. Monitor Component

**Purpose**: Real-time monitoring dashboard and progress reporting

**Module**: `src/pyecod_mini/production/monitor.py`

**Key Features**:
- Live terminal dashboard (using Rich library)
- Progress bars
- Success/failure rates
- Estimated completion time
- Quality metrics distribution

**Implementation**:

```python
from rich.console import Console
from rich.table import Table
from rich.live import Live
from rich.progress import Progress, TaskID

@dataclass
class ProcessingStats:
    """Real-time processing statistics"""
    total_proteins: int
    completed: int
    failed: int
    pending: int
    running: int

    # Quality metrics
    avg_coverage: float = 0.0
    avg_domains: float = 0.0
    quality_distribution: Dict[str, int] = field(default_factory=dict)

    # Performance
    avg_time_per_protein: float = 0.0
    estimated_completion: Optional[datetime] = None

    # Updates
    last_update: datetime = field(default_factory=datetime.now)

class ProductionMonitor:
    """Real-time monitoring of production processing"""

    def __init__(self, tracking_db: TrackingDatabase):
        self.tracking_db = tracking_db
        self.console = Console()

    def get_stats(self) -> ProcessingStats:
        """Get current processing statistics from tracking DB"""
        ...

    def watch(self, refresh_interval: int = 10) -> None:
        """
        Live monitoring dashboard

        Displays:
        - Progress bars
        - Success/failure counts
        - Quality metrics
        - Performance stats
        - Recent completions
        """
        with Live(self._generate_dashboard(), refresh_per_second=1/refresh_interval) as live:
            while True:
                time.sleep(refresh_interval)
                live.update(self._generate_dashboard())

    def _generate_dashboard(self) -> Table:
        """Generate Rich table for dashboard"""
        stats = self.get_stats()

        table = Table(title="pyECOD Mini Production Monitor")

        # Add rows for different metrics
        table.add_row("Total Proteins", str(stats.total_proteins))
        table.add_row("Completed", f"{stats.completed} ({stats.completed/stats.total_proteins*100:.1f}%)")
        table.add_row("Failed", str(stats.failed))
        table.add_row("Running", str(stats.running))
        table.add_row("Pending", str(stats.pending))

        # Quality metrics
        table.add_row("Avg Coverage", f"{stats.avg_coverage:.1%}")
        table.add_row("Avg Domains", f"{stats.avg_domains:.1f}")

        # Performance
        if stats.estimated_completion:
            table.add_row("Est. Completion", stats.estimated_completion.strftime("%Y-%m-%d %H:%M"))

        return table

    def generate_report(self, output_path: Path) -> None:
        """Generate comprehensive processing report"""
        ...
```

**Usage**:
```bash
# Live dashboard
pyecod-mini production monitor --watch

# One-time stats
pyecod-mini production monitor --stats

# Generate report
pyecod-mini production monitor --report output.json
```

### 5. Database Importer

**Purpose**: Import validated results to PostgreSQL production database

**Module**: `src/pyecod_mini/production/database.py`

**Key Features**:
- Collision detection
- Quality filtering
- Batch imports
- Transaction safety
- Verification

**Implementation**:

```python
@dataclass
class ImportResult:
    """Result of database import"""
    attempted: int
    imported: int
    skipped: int
    failed: int
    collisions_detected: int

    details: List[str] = field(default_factory=list)

class DatabaseImporter:
    """Import mini results to PostgreSQL production database"""

    def __init__(self, db_config: DatabaseConfig):
        self.db_config = db_config
        self.conn = self._connect()

    def import_batch(
        self,
        results: List[PartitionResult],
        collision_strategy: str = "skip",  # skip, overwrite, version
        quality_filter: bool = True,
        dry_run: bool = False
    ) -> ImportResult:
        """
        Import batch of results with collision handling

        Args:
            results: Partition results to import
            collision_strategy: How to handle existing data
            quality_filter: Only import good quality results
            dry_run: Don't actually import, just report

        Returns:
            ImportResult with statistics
        """
        ...

    def detect_collisions(self, results: List[PartitionResult]) -> List[str]:
        """Detect proteins that already exist in database"""
        ...

    def verify_import(self, protein_ids: List[str]) -> bool:
        """Verify that proteins were imported correctly"""
        ...
```

## Operational Workflows

### Workflow 1: Initial Production Run

```bash
# 1. Scan batches for representative proteins
pyecod-mini production scan \\
    --batch-pattern "ecod_batch_*" \\
    --reps-only \\
    --output tasks.json

# 2. Submit jobs (start with small test)
pyecod-mini production process \\
    --tasks tasks.json \\
    --max-proteins 100 \\
    --max-concurrent 10

# 3. Monitor progress
pyecod-mini production monitor --watch

# 4. Check results
pyecod-mini production stats

# 5. Import to database (dry run first)
pyecod-mini production import \\
    --check-collisions \\
    --quality-filter \\
    --dry-run

# 6. Actual import
pyecod-mini production import \\
    --collision-strategy skip \\
    --quality-filter

# 7. Verify
pyecod-mini production verify --count 100
```

### Workflow 2: Incremental Updates

```bash
# Process only new proteins
pyecod-mini production process \\
    --skip-existing \\
    --reps-only \\
    --max-concurrent 50
```

### Workflow 3: Retry Failed Jobs

```bash
# Find and retry failed jobs
pyecod-mini production retry \\
    --max-retries 2 \\
    --delay 300
```

## Configuration

**config/production.yml**:
```yaml
database:
  host: dione
  port: 45000
  database: ecod_protein
  user: ecod
  password_env: ECOD_DB_PASSWORD

paths:
  batch_base_dir: /data/ecod/pdb_updates/batches
  tracking_db: /data/ecod/mini_v2/tracking.db
  log_dir: /data/ecod/mini_v2/logs

slurm:
  partition: All
  time: "1:00:00"
  memory: 4G
  cpus_per_task: 1
  max_concurrent_jobs: 50
  job_script_dir: /tmp/pyecod_mini_jobs
  env_setup_script: /home/user/.bashrc

processing:
  reps_only: true
  skip_existing: true
  max_retries: 2
  retry_delay: 300

quality:
  min_coverage: 0.5
  min_confidence: 0.6
  require_reference_lengths: true
  import_quality_threshold: "good"  # good, questionable, all
```

## Monitoring & Alerts

### Metrics to Track

1. **Throughput**: Proteins processed per hour
2. **Success Rate**: % jobs completed successfully
3. **Quality Distribution**: good/questionable/poor breakdown
4. **Performance**: Average processing time
5. **Resource Usage**: Memory, CPU utilization

### Alerts

1. **High Failure Rate**: > 10% failures
2. **Slow Progress**: < 100 proteins/hour
3. **Low Quality**: > 30% poor quality results
4. **Resource Issues**: Memory exhaustion, disk full

## Testing Strategy

### Unit Tests
- Scanner logic
- SLURM script generation
- Database queries
- Import collision detection

### Integration Tests
- End-to-end workflow with 10 proteins
- SLURM submission and monitoring
- Database import and verification

### Production Tests
- 100-protein batch
- Concurrent processing
- Retry mechanism
- Quality filtering

## Deployment Plan

1. **Development Environment** (1 day)
   - Set up on test cluster
   - Process 10 proteins
   - Validate all components

2. **Staging Environment** (2 days)
   - Process 1000 proteins
   - Performance tuning
   - Quality validation

3. **Production Deployment** (3 days)
   - Process all ~40k proteins
   - Real-time monitoring
   - Quality assessment
   - Database import

## Success Criteria

- ✅ Process 40,000 proteins in < 48 hours
- ✅ Success rate >= 95%
- ✅ Average processing time < 10 seconds/protein
- ✅ Quality distribution: >= 60% good, <= 10% poor
- ✅ No data corruption or collisions
- ✅ All results importable to database

---

**Note**: This design builds on proven mini algorithm and avoids over-engineering that failed in the original pyECOD.
