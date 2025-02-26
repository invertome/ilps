#!/bin/bash
# Purpose: Generate static HTML site for offline access to results
# Inputs: analysis/features_final.csv, analysis/predictions.csv, analysis/novel_candidates.csv, candidates/*_blast.out, candidates/*_interpro.tsv, clades_*
# Outputs: output/*/*.csv, output/*/*.png, output/*/ilps.fasta, output/figures/*.html, output/index.html, output/type_pages/*.html, output/summary.html
# Config: config.yaml (max_cpus, ngl_viewer_path)
# Log: pipeline.log (progress, errors)

max_cpus=$(yq e '.max_cpus' config.yaml)
ngl_viewer_url="https://unpkg.com/ngl@2.0.0/dist/ngl.js"  # URL for NGL Viewer
start_time=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting 06_generate_outputs.sh" >> pipeline.log

# Skip if already done
if [ -f output/.done ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping 06_generate_outputs.sh (already done)" >> pipeline.log
    exit 0
fi

# Create output directories
mkdir -p output/prepro output/pro output/mature output/figures output/type_pages output/ngl_viewer

# Download NGL Viewer if not present
if [ ! -f "output/ngl_viewer/ngl.js" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Downloading NGL Viewer" >> pipeline.log
    wget "$ngl_viewer_url" -O "output/ngl_viewer/ngl.js" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to download NGL Viewer" >> pipeline.log; exit 1; }
fi

# Generate tables, plots, FASTA, and metadata (assumes these Python scripts exist)
for type in prepro pro mature; do
    python generate_tables.py analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_blast.out candidates/[0-9]*_interpro.tsv clades_ete_${type}/ clades_autophy_${type}/ clades_treecluster_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_tables.py failed for ${type}" >> pipeline.log; exit 1; }
    python generate_plots.py analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv clades_ete_${type}/ clades_autophy_${type}/ clades_treecluster_${type}/ output/${type} || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_plots.py failed for ${type}" >> pipeline.log; exit 1; }
done

python generate_output_fasta_and_metadata.py analysis/ilp_candidates.fasta analysis/ref_ILPs_filtered.fasta analysis/features_final.csv analysis/predictions.csv analysis/novel_candidates.csv candidates/[0-9]*_interpro.tsv input/[0-9]*_*.fasta output/ || { echo "$(date '+%Y-%m-%d %H:%M:%S') - generate_output_fasta_and_metadata.py failed" >> pipeline.log; exit 1; }

# Generate interactive 3D visualizations with NGL Viewer
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating 3D visualizations" >> pipeline.log
for pdb in analysis/pdbs/*.pdb; do
    base=$(basename "$pdb" .pdb)
    html_file="output/figures/${base}.html"
    cat <<EOF > "$html_file"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>${base} Structure</title>
    <link rel="stylesheet" href="../styles.css">
    <script src="../ngl_viewer/ngl.js"></script>
</head>
<body>
    <div id="viewport" style="width:800px; height:600px;"></div>
    <script>
        document.addEventListener("DOMContentLoaded", function () {
            var stage = new NGL.Stage("viewport");
            stage.loadFile("../analysis/pdbs/${base}.pdb", { defaultRepresentation: true });
        });
    </script>
</body>
</html>
EOF
done

# Generate static HTML site with CSS
echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating static HTML site" >> pipeline.log

# Create styles.css
cat <<EOF > output/styles.css
body {
    font-family: Arial, sans-serif;
    margin: 20px;
}
h1, h2 {
    color: #333;
}
img {
    max-width: 100%;
    height: auto;
}
ul {
    list-style: none;
    padding: 0;
}
li {
    margin: 10px 0;
}
a {
    text-decoration: none;
    color: #007bff;
}
a:hover {
    text-decoration: underline;
}
EOF

# Create index.html
cat <<EOF > output/index.html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ILP Pipeline Results Explorer</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <h1>ILP Pipeline Results Explorer</h1>
    <p>Select a sequence type to explore results:</p>
    <ul>
        <li><a href="type_pages/prepro.html">Prepro</a></li>
        <li><a href="type_pages/pro.html">Pro</a></li>
        <li><a href="type_pages/mature.html">Mature</a></li>
    </ul>
    <p><a href="summary.html">View Metadata Summary</a></p>
</body>
</html>
EOF

# Create type-specific pages
for type in prepro pro mature; do
    cat <<EOF > "output/type_pages/${type}.html"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>${type^} Results</title>
    <link rel="stylesheet" href="../styles.css">
</head>
<body>
    <h1>${type^} Results</h1>
    <h2>Tables</h2>
    <ul>
EOF

    for table in $(ls output/${type}/*.csv 2>/dev/null); do
        base=$(basename "$table")
        echo "        <li><a href=\"../${type}/${base}\">${base}</a></li>" >> "output/type_pages/${type}.html"
    done

    echo "    </ul>" >> "output/type_pages/${type}.html"
    echo "    <h2>Plots</h2>" >> "output/type_pages/${type}.html"
    echo "    <ul>" >> "output_type_pages/${type}.html"

    for plot in $(ls output/${type}/*.png 2>/dev/null); do
        base=$(basename "$plot")
        echo "        <li><img src=\"../${type}/${base}\" alt=\"${base}\"></li>" >> "output_type_pages/${type}.html"
    done

    echo "    </ul>" >> "output_type_pages/${type}.html"
    echo "    <h2>Trees</h2>" >> "output_type_pages/${type}.html"
    echo "    <ul>" >> "output_type_pages/${type}.html"

    for tree in $(ls analysis/${type}_*.tre 2>/dev/null); do
        base=$(basename "$tree")
        echo "        <li><a href=\"../analysis/${base}\">${base}</a></li>" >> "output_type_pages/${type}.html"
    done

    echo "    </ul>" >> "output_type_pages/${type}.html"
    echo "    <h2>Alignments</h2>" >> "output_type_pages/${type}.html"
    echo "    <ul>" >> "output_type_pages/${type}.html"

    for alignment in $(ls analysis/aligned_${type}/*.fasta 2>/dev/null); do
        base=$(basename "$alignment")
        echo "        <li><a href=\"../analysis/aligned_${type}/${base}\">${base}</a></li>" >> "output_type_pages/${type}.html"
    done

    echo "    </ul>" >> "output_type_pages/${type}.html"
    echo "    <h2>3D Structures</h2>" >> "output_type_pages/${type}.html"
    echo "    <ul>" >> "output_type_pages/${type}.html"

    for structure in $(ls output/figures/${type}_*.html 2>/dev/null); do
        base=$(basename "$structure")
        echo "        <li><a href=\"../figures/${base}\">${base}</a></li>" >> "output_type_pages/${type}.html"
    done

    echo "    </ul>" >> "output_type_pages/${type}.html"
    echo "    <a href=\"../index.html\">Back to Home</a>" >> "output_type_pages/${type}.html"
    echo "</body>" >> "output_type_pages/${type}.html"
    echo "</html>" >> "output_type_pages/${type}.html"
done

# Create summary.html
cat <<EOF > output/summary.html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Metadata Summary</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <h1>Metadata Summary</h1>
    <p>Number of Candidates: $(wc -l < analysis/all_candidates.fasta)</p>
    <p>Average Prediction Score: $(awk -F',' '{sum+=$2} END {print sum/NR}' analysis/predictions.csv)</p>
    <!-- Add more summary statistics as needed -->
    <a href="index.html">Back to Home</a>
</body>
</html>
EOF

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "$(date '+%Y-%m-%d %H:%M:%S') - 06_generate_outputs.sh completed in ${runtime}s" >> pipeline.log
touch output/.done
