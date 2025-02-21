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