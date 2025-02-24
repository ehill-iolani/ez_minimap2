# ez_minimap2
This is a test version of a React+Flask app that conducts an alignment using minimap2 and indexes the alignment with samtools.
The results of the alignment can be viewed via the integrated IGV viewer within the page.

The alignment can then be downloaded as a .zip.
Contents of the .zip file are as follows:
- alignment.bam (binary sam file)
- alignment.sam (minimap output)
- alignment.sorted.bam (sorted bam)
- alignment.sorted.bam.bai (bam index)
- reference.fa (uploaded reference)
- target.fa (uploaded targets)

# Installation and usage
## Installation via Docker
The web-app has been containerized using docker and can be built and deployed using the following commands:
```
git clone https://github.com/ehill-iolani/ez_minimap2.git
cd ez_minimap2
docker compose up -d
```

Once the app is up and running navigate to http://localhost:3000.

When you are done using the app and want to shut it down use the following command:
```
docker compose down
```

## Installation and usage for npm + python/flask
For those of you who like to make you lives just a little more difficult (or would like to participate in active development), here are the installation process for installation using npm + python/flask. In an attempt to not introduce conflicts with existing packges, I would highly recommend using a virtual env for this either via venv or through conda. I have provided instructions for making and using a conda virtual env for this example:
```
git clone https://github.com/ehill-iolani/ez_minimap2.git
cd ez_minimap2

# Frontend
npm build
npm start

# Backend
# Create virtual env
conda create -n ez_minimap2 python=3.8
conda activate ez_minimap2

# Installing alignment tools
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels pypi
conda install conda-libmamba-solver
conda config --set solver libmamba
conda install -y bioconda::minimap2
conda install -y bioconda::samtools
conda install -y conda-forge::ncurses

# Install flask and CORS
pip install flask flask-cors

# Run the backend
flask run --port=8080

```