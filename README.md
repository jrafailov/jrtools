# jrtools

Random helpers for the mskilab ecosystem - A collection of utility functions for bioinformatics analysis and system monitoring.

## Installation

You can install the development version of jrtools from GitHub:

```r
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install jrtools
devtools::install_github("jrafailov/jrtools")
```

## Functions

### System Monitoring

#### `hunt(huntdown = NULL)`
Assess local memory usage on mskilab servers. This function provides a tree view of running processes and their memory consumption.

**Parameters:**
- `huntdown`: Optional username to filter results for a specific user

**Returns:** 
- If `huntdown` is NULL: A data.table of all users and their total memory usage in GB
- If `huntdown` is specified: A detailed data.table of that user's processes with PID, status, and memory usage

**Examples:**
```r
library(jrtools)

# Get memory usage for all users
hunt()

# Get detailed process information for a specific user
hunt("username")
```

### Data Processing

#### `concat_file_paths(filepaths, cores = 1)`
Parallelize concatenation of data.tables stored as RDS files. Useful for combining large datasets split across multiple files.

**Parameters:**
- `filepaths`: Vector of RDS file paths containing data.tables
- `cores`: Number of cores to use for parallel processing (default: 1)

**Returns:** Combined data.table from all input files

**Examples:**
```r
# Combine multiple RDS files containing data.tables
files <- c("data1.rds", "data2.rds", "data3.rds")
combined_data <- concat_file_paths(files, cores = 4)
```

### VCF Processing

#### `parsesnpeff(vcf, ...)`
Parse VCF files with SnpEff annotations, with special support for methylation analysis fields. This function processes VCF files and extracts mutation annotations, allelic depth information, and methylation-related metadata.

**Key Features:**
- Supports both standard VCF processing and methylation analysis workflows
- Handles multiple genotype formats (AD, AU/GU/CU/TU/TAR/TIR)
- Extracts SnpEff annotations including gene, transcript, and consequence information
- Supports filtering for coding variants only
- Includes special handling for methylation analysis fields (CpG, METH_PROB, V_HAT, etc.)

### Utility Functions

#### `rand.string(n = 1, length = 12)`
Generate random alphanumeric strings.

**Parameters:**
- `n`: Number of strings to generate (default: 1)
- `length`: Length of each string (default: 12)

**Returns:** Vector of random strings

#### `try2(expr, ..., finally)`
Robust wrapper around `tryCatch` that works well with parallel processing functions.

#### `normalize_path(path)`
Utility function for path normalization that handles special cases like `/dev/null`.

## Dependencies

- `data.table`: For efficient data manipulation
- `ps`: For process monitoring
- `magrittr`: For pipe operator support

Additional suggested packages for full functionality:
- `VariantAnnotation`: For VCF processing
- `parallel`: For multicore processing
- `stringr`: For string operations

## Additional Tools

### SLURM Management Scripts (`scripts/slurmtools`)

A comprehensive collection of bash functions and aliases for managing SLURM jobs on HPC clusters. These tools enhance the standard SLURM commands with additional functionality for monitoring and managing computational workflows.

**Key Features:**
- **Queue Monitoring**: Enhanced `squeue` aliases with custom formatting and real-time watching
  - `sqwatch`: Watch all user jobs with detailed formatting
  - `sqrun`: Monitor only running jobs
  - `sqpend`: Monitor pending jobs
  - `squeuelab`: Display queue status for all lab members
- **Job Management**: 
  - `scancelgrep`: Cancel jobs matching a pattern with confirmation
  - `reassign_jobs`: Reassign jobs to different partitions (useful for Nextflow workflows)
- **Job Inspection**:
  - `job_stderr`: Follow stderr output of running jobs in real-time
  - `job_cwd`: Change to the working directory of a specific job
- **Resource Monitoring**:
  - Enhanced `hunt` function for memory and CPU usage analysis
  - Per-user resource summaries with core-equivalent calculations

**Usage:**
```bash
# Source the tools in your bash profile
source /path/to/jrtools/scripts/slurmtools

# Monitor your running jobs
sqrun

# Cancel all jobs matching a pattern
scancelgrep "my_analysis"

# Follow stderr of a specific job
job_stderr 12345

# Check memory usage across the cluster
hunt --mem
```

### VS Code Configuration (`vscode-config/`)

Curated VS Code settings and extensions optimized for bioinformatics and data science workflows, particularly on remote HPC systems.

**Configuration Files:**
- `settings.json`: Core editor settings with R/Python optimization
- `keybindings.json`: Custom keyboard shortcuts for efficient coding
- `vscode-extensions.json`: Recommended extensions list for bioinformatics
- `nyu___settings.json`: NYU-specific configuration settings

**Key Extensions Included:**
- **Language Support**: R, Python, MATLAB, LaTeX, Groovy
- **Data Science**: CSV editing, PDF viewing, debugging tools
- **Development**: GitHub integration, Copilot, GitLens
- **Remote Work**: SSH remote development tools
- **Themes**: Multiple dark themes optimized for long coding sessions

**Setup:**
```bash
# Copy settings to your VS Code config
cp vscode-config/settings.json ~/.config/Code/User/
cp vscode-config/keybindings.json ~/.config/Code/User/

# Install recommended extensions (requires VS Code CLI)
code --install-extension-list vscode-config/vscode-extensions.json
```

## Use Cases

This package is particularly useful for:

1. **System Administration**: Monitor memory usage and processes on shared computing servers
2. **Bioinformatics Pipelines**: Process and combine large genomic datasets
3. **VCF Analysis**: Parse mutation calls with rich annotation data
4. **Methylation Analysis**: Special support for TAPS and other methylation sequencing workflows
5. **HPC Cluster Management**: Streamlined SLURM job monitoring and management
6. **Remote Development**: Optimized VS Code setup for bioinformatics workflows

## License

MIT License

## Author

Johnathan Rafailov (jrafailov@nygenome.org)

## Contributing

This package is designed for internal use within the mskilab ecosystem, but contributions and suggestions are welcome. Please feel free to open issues or submit pull requests.