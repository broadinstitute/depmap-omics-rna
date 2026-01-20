#!/bin/zsh

set -euo pipefail

uv export \
  --format requirements.txt \
  --all-extras \
  --no-dev \
  --no-hashes \
  --no-editable \
  --no-emit-project \
  > requirements.txt

echo "Updating GCP Function"
gcloud functions deploy depmap-omics-rna \
  --gen2 \
  --runtime="python313" \
  --region="us-central1" \
  --source=. \
  --run-service-account="omics-pipeline-runner@depmap-omics.iam.gserviceaccount.com" \
  --entry-point="run" \
  --trigger-topic="run-depmap-omics-rna" \
  --timeout=540 \
  --memory="4GB" \
  --cpu=1 \
  --max-instances=1 \
  --concurrency=1
