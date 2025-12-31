#!/bin/bash
# ==============================================================================
# Test Script for Bioinformatics R 4.5.1 Container
# ==============================================================================
#
# This script runs comprehensive tests on the container to verify:
#   - Container integrity
#   - R installation
#   - Package availability
#   - Basic functionality
#
# Usage:
#   bash scripts/test_container.sh [CONTAINER_PATH]
#
# If CONTAINER_PATH is not provided, script will look for SIF in project root
#
# ==============================================================================

set -e  # Exit on error

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Default container path
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DEFAULT_SIF="${PROJECT_ROOT}/bioinformatic_r_4.5.1.sif"

# Container to test
CONTAINER="${1:-$DEFAULT_SIF}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

print_header() {
    echo ""
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}======================================================================${NC}"
}

print_test() {
    echo -e "${YELLOW}[TEST $TESTS_RUN] $1${NC}"
}

print_pass() {
    echo -e "${GREEN}  ✓ PASS: $1${NC}"
    ((TESTS_PASSED++))
}

print_fail() {
    echo -e "${RED}  ✗ FAIL: $1${NC}"
    ((TESTS_FAILED++))
}

print_info() {
    echo -e "  ℹ INFO: $1"
}

run_test() {
    ((TESTS_RUN++))
    print_test "$1"
}

# Determine which Singularity command to use
get_singularity_cmd() {
    if command -v singularity &> /dev/null; then
        echo "singularity"
    elif command -v apptainer &> /dev/null; then
        echo "apptainer"
    else
        echo ""
    fi
}

# ------------------------------------------------------------------------------
# PRE-FLIGHT CHECKS
# ------------------------------------------------------------------------------

print_header "Container Test Suite"

echo "Container: $CONTAINER"
echo ""

# Check if Singularity/Apptainer is available
SING_CMD=$(get_singularity_cmd)
if [ -z "$SING_CMD" ]; then
    echo -e "${RED}ERROR: Neither Singularity nor Apptainer found${NC}"
    echo "Please install Singularity or Apptainer to run tests"
    exit 1
fi

echo "Using command: $SING_CMD"
echo "Version: $($SING_CMD --version)"
echo ""

# Check if container exists
if [ ! -f "$CONTAINER" ]; then
    echo -e "${RED}ERROR: Container not found at: $CONTAINER${NC}"
    echo ""
    echo "Usage: bash scripts/test_container.sh [CONTAINER_PATH]"
    exit 1
fi

print_info "Container size: $(du -h "$CONTAINER" | cut -f1)"
echo ""

# ------------------------------------------------------------------------------
# TEST 1: Container Integrity
# ------------------------------------------------------------------------------

run_test "Container Integrity"

if $SING_CMD inspect "$CONTAINER" &> /dev/null; then
    print_pass "Container file is valid"
else
    print_fail "Container file is corrupted or invalid"
fi

# ------------------------------------------------------------------------------
# TEST 2: R Installation
# ------------------------------------------------------------------------------

run_test "R Installation"

R_VERSION=$($SING_CMD exec "$CONTAINER" R --version 2>&1 | head -n1)
if echo "$R_VERSION" | grep -q "R version 4.5.1"; then
    print_pass "R version 4.5.1 is installed"
    print_info "$R_VERSION"
else
    print_fail "Expected R version 4.5.1, got: $R_VERSION"
fi

# ------------------------------------------------------------------------------
# TEST 3: Package Count
# ------------------------------------------------------------------------------

run_test "Package Installation Count"

PACKAGE_COUNT=$($SING_CMD exec "$CONTAINER" \
    Rscript -e "cat(length(.packages(all.available=TRUE)))" 2>/dev/null)

if [ "$PACKAGE_COUNT" -ge 550 ]; then
    print_pass "Found $PACKAGE_COUNT packages (expected ~560)"
else
    print_fail "Found only $PACKAGE_COUNT packages (expected ~560)"
fi

# ------------------------------------------------------------------------------
# TEST 4: Critical Packages
# ------------------------------------------------------------------------------

run_test "Critical Package Availability"

# List of critical packages to test
CRITICAL_PACKAGES=(
    "Seurat"
    "DESeq2"
    "DiffBind"
    "ChIPseeker"
    "monocle3"
    "SingleR"
    "GenomicRanges"
    "ggplot2"
    "dplyr"
    "Biostrings"
)

MISSING_PACKAGES=()

for pkg in "${CRITICAL_PACKAGES[@]}"; do
    if $SING_CMD exec "$CONTAINER" \
        Rscript -e "if (!requireNamespace('$pkg', quietly=TRUE)) quit(status=1)" 2>/dev/null; then
        print_info "✓ $pkg"
    else
        MISSING_PACKAGES+=("$pkg")
        print_info "✗ $pkg (MISSING)"
    fi
done

if [ ${#MISSING_PACKAGES[@]} -eq 0 ]; then
    print_pass "All critical packages available"
else
    print_fail "Missing packages: ${MISSING_PACKAGES[*]}"
fi

# ------------------------------------------------------------------------------
# TEST 5: DESeq2 Functionality
# ------------------------------------------------------------------------------

run_test "DESeq2 Basic Functionality"

TEST_SCRIPT=$(cat << 'EOF'
suppressPackageStartupMessages({
    library(DESeq2)
})

# Create minimal test data
counts <- matrix(rpois(100, lambda=10), ncol=4)
colData <- data.frame(condition = factor(rep(c("A","B"), each=2)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~condition)

# Test should complete without error
cat("SUCCESS")
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "DESeq2 can create objects and run basic operations"
else
    print_fail "DESeq2 functionality test failed"
    print_info "Error: $RESULT"
fi

# ------------------------------------------------------------------------------
# TEST 6: Seurat Functionality
# ------------------------------------------------------------------------------

run_test "Seurat Basic Functionality"

TEST_SCRIPT=$(cat << 'EOF'
suppressPackageStartupMessages({
    library(Seurat)
})

# Create minimal test data
counts <- matrix(rpois(1000, lambda=5), nrow=100)

# Create Seurat object
obj <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

# Test should complete without error
cat("SUCCESS")
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "Seurat can create objects and run basic operations"
else
    print_fail "Seurat functionality test failed"
    print_info "Error: $RESULT"
fi

# ------------------------------------------------------------------------------
# TEST 7: GenomicRanges Functionality
# ------------------------------------------------------------------------------

run_test "GenomicRanges Basic Functionality"

TEST_SCRIPT=$(cat << 'EOF'
suppressPackageStartupMessages({
    library(GenomicRanges)
})

# Create GRanges object
gr <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(100, 200), width = 50)
)

# Test operations
overlaps <- findOverlaps(gr, gr)

cat("SUCCESS")
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "GenomicRanges can create objects and perform operations"
else
    print_fail "GenomicRanges functionality test failed"
    print_info "Error: $RESULT"
fi

# ------------------------------------------------------------------------------
# TEST 8: ggplot2 Plotting
# ------------------------------------------------------------------------------

run_test "ggplot2 Plotting Capability"

TEST_SCRIPT=$(cat << 'EOF'
suppressPackageStartupMessages({
    library(ggplot2)
})

# Create simple plot
data <- data.frame(x = 1:10, y = 1:10)
p <- ggplot(data, aes(x, y)) + geom_point()

# Save to temp file (don't display)
temp_file <- tempfile(fileext = ".pdf")
ggsave(temp_file, p, width = 5, height = 5, device = "pdf")

# Check file was created
if (file.exists(temp_file)) {
    cat("SUCCESS")
    unlink(temp_file)
}
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "ggplot2 can create and save plots"
else
    print_fail "ggplot2 plotting test failed"
    print_info "Error: $RESULT"
fi

# ------------------------------------------------------------------------------
# TEST 9: File I/O Operations
# ------------------------------------------------------------------------------

run_test "File I/O Operations"

TEST_SCRIPT=$(cat << 'EOF'
# Test write
temp_file <- tempfile(fileext = ".csv")
data <- data.frame(x = 1:10, y = 11:20)
write.csv(data, temp_file, row.names = FALSE)

# Test read
data_read <- read.csv(temp_file)

# Verify
if (nrow(data_read) == 10 && ncol(data_read) == 2) {
    cat("SUCCESS")
}

unlink(temp_file)
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "File read/write operations work correctly"
else
    print_fail "File I/O test failed"
fi

# ------------------------------------------------------------------------------
# TEST 10: Memory and Performance
# ------------------------------------------------------------------------------

run_test "Memory Allocation Test"

TEST_SCRIPT=$(cat << 'EOF'
# Allocate ~100MB matrix
big_matrix <- matrix(rnorm(1e7), nrow=1000)

# Perform operation
result <- colMeans(big_matrix)

cat("SUCCESS")
EOF
)

RESULT=$($SING_CMD exec "$CONTAINER" Rscript -e "$TEST_SCRIPT" 2>&1)

if echo "$RESULT" | grep -q "SUCCESS"; then
    print_pass "Can allocate and process large matrices"
else
    print_fail "Memory allocation test failed"
fi

# ------------------------------------------------------------------------------
# TEST SUMMARY
# ------------------------------------------------------------------------------

print_header "Test Summary"

echo ""
echo "Total tests run: $TESTS_RUN"
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    echo ""
    echo "The container is ready for use."
    echo ""
    echo "Quick start:"
    echo "  $SING_CMD exec $CONTAINER R"
    echo ""
    echo "Run a script:"
    echo "  $SING_CMD exec $CONTAINER Rscript your_script.R"
    echo ""
    EXIT_CODE=0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    echo ""
    echo "Please review the failed tests above."
    echo "The container may not function correctly."
    echo ""
    EXIT_CODE=1
fi

exit $EXIT_CODE