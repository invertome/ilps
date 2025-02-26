#!/bin/bash
# Purpose: Step 06 - Generate tables, plots, FASTA files, metadata TSV, interactive 3D visualizations, and launch web app
# Inputs: analysis/features_final.csv, analysis/predictions.csv, analysis/novel_candidates.csv, candidates/*_blast.out, candidates/*_interpro.tsv, clades_*
# Outputs: output/*/*.csv, output/*/*.png, output/*/ilps.fasta, output/comparative_metadata.tsv, output/figures/*.html, web app
# Config: config.yaml (max_cpus, ngl_viewer_path, flask_app_port)
# Log: pipeline.log (progress, errors, profiling)
# Notes: Generates interactive 3D visualizations with NGL Viewer and launches Flask web app

max_cpus=$(yq e '.max_cpus' config.yaml)
ngl_viewer_path=$(yq e '.ngl_viewer_path' config.yaml)
flask_app_port=$(yq e '.flask_app_port' config.yaml)
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 06_generate_outputs.sh" >> pipeline.log
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory before: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"

if [ -f output/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 06_generate_outputs.sh (already done)" >> pipeline.log
    exit 0
fi

mkdir -p output/prepro output/pro output/mature output/figures

# Extract final features including structural data
python extract_features.py analysis/all_candidates.fasta analysis/features_final.csv || { echo "$(date '+%Y-%m-%d %H:%M:%S') - extract_features.py failed" >> pipeline.log; exit 1; }

for type in prepro pro mature; do
    python generate_tables.py analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_${type}/ clades_autophy_${type}/ clades_treecluster_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_tables.py failed for ${type}" >> pipeline.log; exit 1; }
    python generate_plots.py analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_${type}/ clades_autophy_${type}/ clades_treecluster_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_plots.py failed for ${type}" >> pipeline.log; exit 1; }
done

python generate_output_fasta_and_metadata.py analysis/ilp_candidates.fasta analysis/ref_ILPs_filtered.fasta analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_interpro.tsv input/[0-9]*_*.fasta output/ || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_output_fasta_and_metadata.py failed" >> pipeline.log; exit 1; }

# Generate interactive 3D visualizations with NGL Viewer
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating interactive 3D visualizations" >> pipeline.log
for pdb in analysis/pdbs/*.pdb; do
    base=$(basename "$pdb" .pdb)
    html_file="output/figures/${base}.html"
    cat <<EOF > "$html_file"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>${base} Structure</title>
    <script src="${ngl_viewer_path}/ngl.js"></script>
</head>
<body>
    <div id="viewport" style="width:800px; height:600px;"></div>
    <script>
        document.addEventListener("DOMContentLoaded", function () {
            var stage = new NGL.Stage("viewport");
            stage.loadFile("${pdb}", { defaultRepresentation: true });
        });
    </script>
</body>
</html>
EOF
done

# Launch Flask web app to explore results
echo "$(date '+%Y-%m-%d %H:%M:%S') - Launching Flask web app on port $flask_app_port" >> pipeline.log
nohup python results_explorer.py --port "$flask_app_port" > flask.log 2>&1 &

end_time=$(date +%s)
runtime=$((end_time - start_time))
python -c "import psutil; print(f'$(date '+%Y-%m-%d %H:%M:%S') - Memory after: {psutil.virtual_memory().percent}%', file=open('pipeline.log', 'a'))"
echo "$(date '+%Y-%m-%d %H:%M:%S') - 06_generate_outputs.sh completed in ${runtime}s" >> pipeline.log
touch output/.done
