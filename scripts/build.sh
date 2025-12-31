#!/bin/bash
# ==============================================================================
# Build Script for Bioinformatics R 4.5.1 Container
# ==============================================================================
#
# This script builds a Docker image and optionally converts it to Singularity SIF
#
# Usage:
#   bash scripts/build.sh [OPTIONS]
#
# Options:
#   --docker-only    Build Docker image only (skip Singularity conversion)
#   --sif-only       Convert existing Docker image to SIF only
#   --no-push        Don't push to DockerHub
#   --help           Show this help message
#
# Requirements:
#   - Docker installed and running
#   - Singularity/Apptainer (for SIF conversion)
#   - DockerHub account configured (for pushing)
#
# ==============================================================================

set -e  # Exit on error

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Docker image details
DOCKER_USER="gynecoloji"
IMAGE_NAME="bioinformatic_r_4_5_1"
IMAGE_TAG="v1"
FULL_IMAGE="${DOCKER_USER}/${IMAGE_NAME}:${IMAGE_TAG}"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DOCKERFILE="${PROJECT_ROOT}/definitions/docker_file"
SIF_OUTPUT="${PROJECT_ROOT}/bioinformatic_r_4.5.1.sif"

# Build options
BUILD_DOCKER=true
BUILD_SIF=true
PUSH_DOCKER=true

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

print_header() {
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}======================================================================${NC}"
}

print_step() {
    echo -e "${GREEN}>>> $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

print_error() {
    echo -e "${RED}ERROR: $1${NC}"
}

show_help() {
    cat << EOF
Build Script for Bioinformatics R 4.5.1 Container

Usage:
    bash scripts/build.sh [OPTIONS]

Options:
    --docker-only    Build Docker image only (skip Singularity conversion)
    --sif-only       Convert existing Docker image to SIF only
    --no-push        Don't push to DockerHub
    --help           Show this help message

Examples:
    # Build everything (Docker + SIF) and push to DockerHub
    bash scripts/build.sh

    # Build Docker image only
    bash scripts/build.sh --docker-only

    # Convert existing Docker image to SIF
    bash scripts/build.sh --sif-only

    # Build but don't push to DockerHub
    bash scripts/build.sh --no-push

EOF
}

check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running"
        exit 1
    fi
    
    print_step "Docker is available: $(docker --version)"
}

check_singularity() {
    if command -v singularity &> /dev/null; then
        print_step "Singularity is available: $(singularity --version)"
        return 0
    elif command -v apptainer &> /dev/null; then
        print_step "Apptainer is available: $(apptainer --version)"
        return 0
    else
        print_warning "Neither Singularity nor Apptainer found"
        return 1
    fi
}

build_docker_image() {
    print_header "Building Docker Image"
    
    print_step "Building image: ${FULL_IMAGE}"
    print_step "Using Dockerfile: ${DOCKERFILE}"
    
    # Build the Docker image
    docker build \
        -f "$DOCKERFILE" \
        -t "$FULL_IMAGE" \
        "$PROJECT_ROOT"
    
    print_step "Docker image built successfully"
    
    # Show image details
    echo ""
    print_step "Image details:"
    docker images "$FULL_IMAGE"
    echo ""
}

push_docker_image() {
    print_header "Pushing to DockerHub"
    
    print_step "Pushing image: ${FULL_IMAGE}"
    
    # Check if logged in to DockerHub
    if ! docker info | grep -q "Username"; then
        print_warning "Not logged in to DockerHub"
        echo "Please run: docker login"
        read -p "Press Enter to continue after logging in, or Ctrl+C to cancel..."
    fi
    
    # Push the image
    docker push "$FULL_IMAGE"
    
    print_step "Image pushed successfully to DockerHub"
    echo ""
    print_step "Image URL: https://hub.docker.com/r/${DOCKER_USER}/${IMAGE_NAME}"
}

build_singularity_image() {
    print_header "Converting to Singularity SIF"
    
    # Determine which command to use
    local SING_CMD=""
    if command -v singularity &> /dev/null; then
        SING_CMD="singularity"
    elif command -v apptainer &> /dev/null; then
        SING_CMD="apptainer"
    else
        print_error "Neither Singularity nor Apptainer found"
        return 1
    fi
    
    print_step "Using command: ${SING_CMD}"
    print_step "Output SIF: ${SIF_OUTPUT}"
    
    # Remove existing SIF if present
    if [ -f "$SIF_OUTPUT" ]; then
        print_warning "Removing existing SIF file"
        rm -f "$SIF_OUTPUT"
    fi
    
    # Build SIF from Docker image
    print_step "Converting Docker image to SIF..."
    print_step "This may take several minutes..."
    
    $SING_CMD build "$SIF_OUTPUT" "docker-daemon://${FULL_IMAGE}"
    
    print_step "SIF file created successfully"
    
    # Show SIF details
    echo ""
    print_step "SIF file details:"
    ls -lh "$SIF_OUTPUT"
    echo ""
    
    # Test the SIF
    print_step "Testing SIF file..."
    $SING_CMD exec "$SIF_OUTPUT" R --version
    echo ""
}

verify_packages() {
    print_header "Verifying Package Installation"
    
    print_step "Checking number of installed packages..."
    
    # Use docker or singularity depending on what we built
    if [ "$BUILD_DOCKER" = true ]; then
        docker run --rm "$FULL_IMAGE" \
            Rscript -e "cat('Total packages:', length(.packages(all.available=TRUE)), '\n')"
    elif [ -f "$SIF_OUTPUT" ]; then
        if command -v singularity &> /dev/null; then
            singularity exec "$SIF_OUTPUT" \
                Rscript -e "cat('Total packages:', length(.packages(all.available=TRUE)), '\n')"
        else
            apptainer exec "$SIF_OUTPUT" \
                Rscript -e "cat('Total packages:', length(.packages(all.available=TRUE)), '\n')"
        fi
    fi
    
    echo ""
}

# ------------------------------------------------------------------------------
# PARSE COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case $1 in
        --docker-only)
            BUILD_SIF=false
            shift
            ;;
        --sif-only)
            BUILD_DOCKER=false
            PUSH_DOCKER=false
            shift
            ;;
        --no-push)
            PUSH_DOCKER=false
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------------------------

print_header "Bioinformatics R 4.5.1 Container Build Script"

echo ""
echo "Configuration:"
echo "  Docker image: ${FULL_IMAGE}"
echo "  Build Docker: ${BUILD_DOCKER}"
echo "  Build SIF: ${BUILD_SIF}"
echo "  Push to DockerHub: ${PUSH_DOCKER}"
echo ""

# Check prerequisites
print_step "Checking prerequisites..."

if [ "$BUILD_DOCKER" = true ]; then
    check_docker
fi

if [ "$BUILD_SIF" = true ]; then
    if ! check_singularity; then
        print_error "Singularity/Apptainer required for SIF build"
        exit 1
    fi
fi

# Verify Dockerfile exists
if [ ! -f "$DOCKERFILE" ]; then
    print_error "Dockerfile not found at: $DOCKERFILE"
    exit 1
fi

echo ""
read -p "Press Enter to continue with build, or Ctrl+C to cancel..."
echo ""

# Build Docker image
if [ "$BUILD_DOCKER" = true ]; then
    build_docker_image
fi

# Push to DockerHub
if [ "$PUSH_DOCKER" = true ] && [ "$BUILD_DOCKER" = true ]; then
    push_docker_image
fi

# Build Singularity SIF
if [ "$BUILD_SIF" = true ]; then
    # If we didn't build Docker, pull from DockerHub first
    if [ "$BUILD_DOCKER" = false ]; then
        print_step "Pulling Docker image from DockerHub..."
        docker pull "$FULL_IMAGE"
    fi
    build_singularity_image
fi

# Verify package installation
verify_packages

# ------------------------------------------------------------------------------
# BUILD COMPLETE
# ------------------------------------------------------------------------------

print_header "Build Complete!"

echo ""
echo "Summary:"
if [ "$BUILD_DOCKER" = true ]; then
    echo -e "  ${GREEN}✓${NC} Docker image built: ${FULL_IMAGE}"
fi
if [ "$PUSH_DOCKER" = true ]; then
    echo -e "  ${GREEN}✓${NC} Pushed to DockerHub"
fi
if [ "$BUILD_SIF" = true ]; then
    echo -e "  ${GREEN}✓${NC} SIF file created: ${SIF_OUTPUT}"
fi
echo ""

if [ -f "$SIF_OUTPUT" ]; then
    echo "To test the SIF file:"
    echo "  singularity exec ${SIF_OUTPUT} R --version"
    echo ""
    echo "To run an R script:"
    echo "  singularity exec ${SIF_OUTPUT} Rscript your_script.R"
    echo ""
fi

if [ "$PUSH_DOCKER" = true ]; then
    echo "To pull from DockerHub:"
    echo "  singularity pull docker://${FULL_IMAGE}"
    echo ""
fi

echo "See docs/INSTALLATION.md for usage instructions."
echo ""

exit 0