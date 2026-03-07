# IgDiscover.jl Makefile

IMAGE_NAME := igdiscover-parity
LIMIT ?= 200
ITERATIONS ?= 1
RESULTS ?= $(PWD)/results

.PHONY: help build test test-reads test-interactive clean

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

build: ## Build the Docker parity test image
	docker build -t $(IMAGE_NAME) .

test: build ## Run parity test with synthetic reads (LIMIT=200). Results in $(RESULTS)
	@mkdir -p $(RESULTS)
	docker run --rm \
		-v $(RESULTS):/results \
		-e LIMIT=$(LIMIT) \
		-e ITERATIONS=$(ITERATIONS) \
		$(IMAGE_NAME)

test-reads: build ## Run parity test with custom reads: make test-reads READS=/path/to/reads.fasta.gz
	@test -n "$(READS)" || (echo "Usage: make test-reads READS=/path/to/reads.fasta.gz"; exit 1)
	@mkdir -p $(RESULTS)
	docker run --rm \
		-v $(RESULTS):/results \
		-v $(READS):/data/reads.fasta.gz:ro \
		-e LIMIT=$(LIMIT) \
		-e ITERATIONS=$(ITERATIONS) \
		$(IMAGE_NAME)

test-interactive: build ## Start interactive shell in the parity test container
	@mkdir -p $(RESULTS)
	docker run --rm -it \
		-v $(RESULTS):/results \
		$(IMAGE_NAME) bash

test-unit: ## Run Julia unit tests (requires Julia locally)
	julia --project=. -e 'using Pkg; Pkg.test()'

clean: ## Remove Docker image and local results
	docker rmi $(IMAGE_NAME) 2>/dev/null || true
	rm -rf $(RESULTS)/
