import React, { useState, useEffect } from 'react';
import './App.css';
import igv from '../node_modules/igv/dist/igv.esm.min.js';

function App() {
  const [referenceFile, setReferenceFile] = useState(null);
  const [targetFile, setTargetFile] = useState(null);
  const [result, setResult] = useState('');
  const [error, setError] = useState('');
  const [downloadUrl, setDownloadUrl] = useState('');
  const [igvBrowser, setIgvBrowser] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [streamedOutput, setStreamedOutput] = useState('');
  const [isMinION, setIsMinION] = useState(false);
  const [isDownloading, setIsDownloading] = useState(false);

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

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setResult('');
    setDownloadUrl('');
    setStreamedOutput('');
    setIsLoading(true);

    const formData = new FormData();
    formData.append('reference', referenceFile);
    formData.append('target', targetFile);
    formData.append('isMinION', isMinION); // Add the checkbox state to the form data

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

      while (!done) {
        const { value, done: doneReading } = await reader.read();
        done = doneReading;
        const chunk = decoder.decode(value, { stream: true });
        setStreamedOutput((prev) => prev + chunk);
      }

      setIsLoading(false);
      setDownloadUrl('http://localhost:8080/compress_and_download'); // Set the download URL
      loadIgvBrowser();
    } catch (err) {
      setError(err.message);
      setIsLoading(false);
    }
  };

  const handleDownload = async () => {
    setIsDownloading(true);
    try {
      const response = await fetch('http://localhost:8080/compress_and_download', {
        method: 'GET',
      });

      if (!response.ok) {
        throw new Error('Failed to compress and download results');
      }

      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'output_alignment_dl.zip';
      document.body.appendChild(a);
      a.click();
      a.remove();
    } catch (err) {
      setError(err.message);
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
    <div className="App">
      <header className="App-header">
        <h1>EZ Minimap2 Alignment</h1>
        <form onSubmit={handleSubmit}>
          <div>
            <label>
              Reference File:
              <input
                type="file"
                onChange={(e) => setReferenceFile(e.target.files[0])}
                required
              />
            </label>
          </div>
          <div>
            <label>
              Target File:
              <input
                type="file"
                onChange={(e) => setTargetFile(e.target.files[0])}
                required
              />
            </label>
          </div>
          <div>
            <label>
              <input
                type="checkbox"
                checked={isMinION}
                onChange={(e) => setIsMinION(e.target.checked)}
              />
              <span>Are you uploading FASTQs produced by a MinION sequencer?</span>
            </label>
          </div>
          <button type="submit">Run Alignment</button>
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
        {downloadUrl && (
          <button onClick={handleDownload}>Download Results</button>
        )}
        <p className="alignment-logs">Alignment Logs:</p>
        <textarea
          value={streamedOutput}
          readOnly
          style={{ width: '80%', height: '300px', marginTop: '10px' }}
        />
        <div id="igv-div" style={{ width: '80%', height: '1200px' }}></div>
      </header>
    </div>
  );
}

export default App;