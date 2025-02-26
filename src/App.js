import React, { useState, useEffect } from 'react';
import './App.css';
import igv from '../node_modules/igv/dist/igv.esm.min.js';
//import JSZip from 'jszip';

// IMPORTANT: when deploying a server you should use the server's IP address instead of localhost!!

function App() {
  const [referenceFile, setReferenceFile] = useState(null);
  const [targetFile, setTargetFile] = useState(null);
  const [result, setResult] = useState('');
  const [error, setError] = useState('');
  //const [downloadUrl, setDownloadUrl] = useState('');
  const [igvBrowser, setIgvBrowser] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [streamedOutput, setStreamedOutput] = useState('');
  const [isMinION, setIsMinION] = useState(false);
  const [isDownloading, setIsDownloading] = useState(false);
  const [hasCachedOutput, setHasCachedOutput] = useState(false);

  useEffect(() => {
    const deleteOutputs = async () => {
      try {
        console.log('Attempting to delete outputs');
        const response = await fetch('http://localhost:8080/delete_outputs', {
          method: 'DELETE',
        });
        if (!response.ok) {
          throw new Error('Failed to delete outputs');
        }
        console.log('Outputs deleted successfully');
      } catch (err) {
        console.error('Failed to delete outputs:', err);
      }
    };

    deleteOutputs();
  }, []);

  useEffect(() => {
    const cachedOutput = localStorage.getItem('alignmentOutput');
    if (cachedOutput) {
      setHasCachedOutput(true);
    }
  }, []);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setResult('');
    setStreamedOutput('');
    setIsLoading(true);

    const formData = new FormData();
    formData.append('reference', referenceFile);
    formData.append('target', targetFile);
    formData.append('isMinION', isMinION);

    try {
      const response = await fetch('http://localhost:8080/run_alignment', {
        method: 'POST',
        body: formData,
      });

      if (!response.ok) {
        throw new Error('Network response was not ok');
      }

      const reader = response.body.getReader();
      const decoder = new TextDecoder('utf-8');
      let done = false;
      let output = '';

      while (!done) {
        const { value, done: doneReading } = await reader.read();
        done = doneReading;
        const chunk = decoder.decode(value, { stream: true });
        output += chunk;
        setStreamedOutput((prev) => prev + chunk);
      }

      // Store the output in local storage
      localStorage.setItem('alignmentOutput', output);
      setHasCachedOutput(true);

      setIsLoading(false);
      loadIgvBrowser();
    } catch (err) {
      setError(err.message);
      setIsLoading(false);
    }
  };

  const handleDownload = async () => {
    setIsDownloading(true);
    try {
      const response = await fetch('http://localhost:8080/get_output_files', {
        method: 'GET',
      });

      if (!response.ok) {
        throw new Error('Failed to fetch output files');
      }

      const files = await response.json();

      for (const file of files) {
        const fileResponse = await fetch(`http://localhost:8080/output_alignment/${file}`);
        if (!fileResponse.ok) {
          throw new Error(`Failed to fetch file: ${file}`);
        }

        const blob = await fileResponse.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = file;
        document.body.appendChild(a);
        a.click();
        a.remove();
      }
    } catch (err) {
      setError(err.message);
      console.log(err);
    } finally {
      setIsDownloading(false);
    }
  };

  const loadIgvBrowser = async () => {
    const igvDiv = document.getElementById('igv-div');
    const igvOptions = {
      reference: {
        id: 'Output',
        fastaURL: 'http://localhost:8080/output_alignment/reference.fa',
        indexed: false,
        tracks: [
          {
            name: 'Mapping',
            type: 'alignment',
            url: 'http://localhost:8080/output_alignment/alignment.sorted.bam',
            indexURL: 'http://localhost:8080/output_alignment/alignment.sorted.bam.bai',
            displayMode: 'EXPANDED',
            format: 'bam',
            indexed: true,
          },
        ],
      },
    };

    if (igvBrowser) {
      igvBrowser.dispose();
    }

    const browser = await igv.createBrowser(igvDiv, igvOptions);
    setIgvBrowser(browser);
  };

  useEffect(() => {
    // Cleanup IGV browser on component unmount
    return () => {
      if (igvBrowser) {
        igvBrowser.dispose();
      }
    };
  }, [igvBrowser]);

  return (
    <div className="App bg-dark">
      <header className="App-header bg-white d-flex flex-column align-items-center justify-content-center">
        <h1 className="display-4">EZ Minimap2 Alignment</h1>
        <form className="form-group" onSubmit={handleSubmit}>
          <div className="d-flex justify-content-center align-items-center">
            <label className="form-label d-flex align-items-center">
              Reference File:
              <input
                type="file"
                className="form-control"
                onChange={(e) => setReferenceFile(e.target.files[0])}
                required
                accept=".fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz"
              />
            </label>
          </div>
          <div className="d-flex justify-content-center align-items-center">
            <label className="form-label">
              Target File:
              <input
                type="file"
                className="form-control"
                onChange={(e) => setTargetFile(e.target.files[0])}
                required
                accept=".fa,.fasta,.fna,.fastq,.fq,application/gzip,.fa.gz,.fasta.gz,.fna.gz,.fastq.gz,.fq.gz"
              />
            </label>
          </div>
          <div>
            <label className="form-label">
              <input
                type="checkbox"
                className="form-check-input"
                checked={isMinION}
                onChange={(e) => setIsMinION(e.target.checked)}
              />
              <span>Are you uploading FASTQs produced by a MinION sequencer?</span>
            </label>
          </div>
          <div className="d-flex justify-content-center">
            <button type="submit" className="btn btn-primary mt-3">Run Alignment</button>
          </div>
        </form>
        {isLoading && (
          <div className="loader-container">
            <div className="loader"></div>
            <div className="loader-text">Loading...</div>
          </div>
        )}
        {isDownloading && (
          <div className="downloader-container">
            <div className="downloader"></div>
            <div className="downloader-text">Compressing and preparing download...</div>
          </div>
        )}
        {error && <p className="error">{error}</p>}
        {result && <p className="result">{result}</p>}
        {hasCachedOutput && (
          <button onClick={handleDownload} className="btn btn-success mt-3">Download Results</button>
        )}
        <p className="alignment-logs">Alignment Logs:</p>
        <textarea
          value={streamedOutput}
          readOnly
          className="form-control"
          style={{ width: '80%', height: '300px', marginTop: '10px' }}
        />
        <div id="igv-div" className="mt-3" style={{ width: '80%', height: '1200px' }}></div>
      </header>
    </div>
  );
}

export default App;