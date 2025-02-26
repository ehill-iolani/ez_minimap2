import os
import subprocess
import logging
from flask import Flask, request, jsonify, send_file, Response
from flask_cors import CORS
import shutil
import gzip
import time
import uuid
import zipfile

# Create flask app
app = Flask(__name__)

# CORS allows for cross-origin requests (ie requests from a different server)
CORS(app, resources={r"/*": {"origins": "*"}})

# Configure logging
logging.basicConfig(level=logging.DEBUG)

####################################
### Setting the guts of the shit ###
####################################
# The guts alignment minimap2 + samtools
def run_minimap2(reference, target, is_minion, session_id):
    """
    Run minimap2 alignment with the given reference and target sequences.
    
    :param reference: Path to the reference sequence file
    :param target: Path to the target sequence file
    :param is_minion: Boolean indicating if the sequences are from a MinION sequencer
    :param session_id: Unique session ID
    :return: Generator yielding alignment result as a string
    """
    try:
        yield "Starting alignment and processing...\n"
        
        # Create a unique directory for the session
        session_dir = os.path.join('/tmp', session_id)
        os.makedirs(session_dir, exist_ok=True)
        
        # Move reference and target files to the session directory
        reference_path = os.path.join(session_dir, os.path.basename(reference))
        target_path = os.path.join(session_dir, os.path.basename(target))
        shutil.move(reference, reference_path)
        shutil.move(target, target_path)
        
        # Check if the sequences are .gz files and unzip if necessary
        start_time = time.time()
        if reference_path.endswith('.gz'):
            reference_unzipped = reference_path.replace('.gz', '')
            with open(reference_unzipped, 'wb') as f_out, gzip.open(reference_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)
            reference_path = reference_unzipped
            end_time = time.time()
            yield f"Time taken to unzip reference: {end_time - start_time} seconds\n"
            yield "Reference file unzipped successfully\n"
        else:
            yield "Uploaded reference file does not need to be unzipped\n"
        if target_path.endswith('.gz'):
            target_unzipped = target_path.replace('.gz', '')
            with open(target_unzipped, 'wb') as f_out, gzip.open(target_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)
            target_path = target_unzipped
            end_time = time.time()
            yield f"Time taken to unzip targets: {end_time - start_time} seconds\n"
            yield "Target file unzipped successfully\n"
        else:
            yield "Uploaded target files(s) do not need to be unzipped\n"
        
        # Write reference and target sequences to files
        with open(os.path.join(session_dir, 'reference.fa'), 'w') as f:
            with open(reference_path, 'r') as ref:
                f.write(ref.read())
        with open(os.path.join(session_dir, 'target.fa'), 'w') as f:
            with open(target_path, 'r') as tar:
                f.write(tar.read())

        # Run minimap2
        start_time = time.time()
        if is_minion:
            minimap2_cmd = ['minimap2', '-ax', 'map-ont', '-t', '8', reference_path, target_path, '-o', os.path.join(session_dir, 'alignment.sam')]
        else:
            minimap2_cmd = ['minimap2', '-asm5', '-t', '8', reference_path, target_path, '-o', os.path.join(session_dir, 'alignment.sam')]
        with subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                yield line
            for line in proc.stderr:
                yield line
            proc.wait()
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, minimap2_cmd)
        end_time = time.time()
        yield f"Time taken to run minimap2: {end_time - start_time} seconds\n"
        yield "minimap2 command ran successfully\n"
        
        # Convert SAM to BAM
        start_time = time.time()
        samtools_view_cmd = ['samtools', 'view', '-@ 8', '-bS', os.path.join(session_dir, 'alignment.sam'), '-o', os.path.join(session_dir, 'alignment.bam')]
        with subprocess.Popen(samtools_view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                yield line
            for line in proc.stderr:
                yield line
            proc.wait()
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, samtools_view_cmd)
        end_time = time.time()
        yield f"Time taken to convert SAM to BAM: {end_time - start_time} seconds\n"
        yield "samtools view command ran successfully\n"
        
        # Sort BAM file
        start_time = time.time()
        samtools_sort_cmd = ['samtools', 'sort', '-@ 8', os.path.join(session_dir, 'alignment.bam'), '-o', os.path.join(session_dir, 'alignment.sorted.bam')]
        with subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                yield line
            for line in proc.stderr:
                yield line
            proc.wait()
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, samtools_sort_cmd)
        end_time = time.time()
        yield f"Time taken to sort BAM file: {end_time - start_time} seconds\n"
        yield "samtools sort command ran successfully\n"
        
        # Index BAM file
        start_time = time.time()
        samtools_index_cmd = ['samtools', 'index', '-@ 8', os.path.join(session_dir, 'alignment.sorted.bam')]
        with subprocess.Popen(samtools_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                yield line
            for line in proc.stderr:
                yield line
            proc.wait()
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, samtools_index_cmd)
        end_time = time.time()
        yield f"Time taken to index BAM file: {end_time - start_time} seconds\n"
        yield "samtools index command ran successfully\n"
        
        # Move results to session directory
        yield "Results are being moved...please be patient\n"
        start_time = time.time()
        output_dir = os.path.join(session_dir, 'output_alignment')
        os.makedirs(output_dir, exist_ok=True)
        shutil.move(os.path.join(session_dir, 'alignment.sorted.bam'), os.path.join(output_dir, 'alignment.sorted.bam'))
        shutil.move(os.path.join(session_dir, 'alignment.sorted.bam.bai'), os.path.join(output_dir, 'alignment.sorted.bam.bai'))
        shutil.move(os.path.join(session_dir, 'alignment.sam'), os.path.join(output_dir, 'alignment.sam'))
        shutil.move(os.path.join(session_dir, 'alignment.bam'), os.path.join(output_dir, 'alignment.bam'))
        shutil.move(os.path.join(session_dir, 'reference.fa'), os.path.join(output_dir, 'reference.fa'))
        shutil.move(os.path.join(session_dir, 'target.fa'), os.path.join(output_dir, 'target.fa'))
        end_time = time.time()
        yield f"Time taken to move results: {end_time - start_time} seconds\n"
        yield "Results moved successfully\n"
        
        yield "Alignment and processing completed successfully!\nResults are ready to be compressed.\n"

    except subprocess.CalledProcessError as e:
        yield f"An error occurred: {e.stderr}\n"

# Compressing the results
def compress_results(session_id):
    try:
        start_time = time.time()
        session_dir = os.path.join('/tmp', session_id)
        output_dir = os.path.join(session_dir, 'output_alignment')
        with zipfile.ZipFile(os.path.join(session_dir, 'output_alignment_dl.zip'), 'w', compression=zipfile.ZIP_STORED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), output_dir))
        end_time = time.time()
        logging.debug(f"Time taken to compress results: {end_time - start_time} seconds")
    except Exception as e:
        logging.error(f"Error compressing results: {e}")
        raise

##############################
### Define routes and shit ###
##############################
# Route to upload and run the alignment
@app.route('/run_alignment', methods=['POST'])
def run_alignment():
    app.logger.debug('Received request for run_alignment')
    if 'reference' not in request.files or 'target' not in request.files:
        app.logger.error('Reference and target files are required')
        return jsonify({'error': 'Reference and target files are required'}), 400

    session_id = request.form.get('sessionId')
    if not session_id:
        session_id = str(uuid.uuid4())

    reference = request.files['reference']
    target = request.files['target']
    is_minion = request.form.get('isMinION') == 'true'

    reference_path = os.path.join('/tmp', reference.filename)
    target_path = os.path.join('/tmp', target.filename)

    reference.save(reference_path)
    target.save(target_path)

    app.logger.debug(f'Saved reference file to {reference_path}')
    app.logger.debug(f'Saved target file to {target_path}')

    return Response(run_minimap2(reference_path, target_path, is_minion, session_id), content_type='text/plain')

# Route to intiate the file compression and download
@app.route('/compress_and_download', methods=['GET'])
def compress_and_download():
    session_id = request.args.get('sessionId')
    if not session_id:
        return jsonify({'error': 'Session ID is required'}), 400

    try:
        compress_results(session_id)
        session_dir = os.path.join('/tmp', session_id)
        return send_file(os.path.join(session_dir, 'output_alignment_dl.zip'), as_attachment=True)
    except Exception as e:
        app.logger.error(f"Error compressing or downloading results: {e}")
        return jsonify({'error': 'Error compressing or downloading results'}), 500

# Route to visualize the alignment
@app.route('/output_alignment/<session_id>/<path:filename>', methods=['GET'])
def download_file(session_id, filename):
    try:
        file_path = os.path.join('/tmp', session_id, 'output_alignment', filename)
        if os.path.exists(file_path):
            return send_file(file_path, as_attachment=True)
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        app.logger.error(f"File not found: {filename}")
        return jsonify({'error': 'File not found'}), 404
    except Exception as e:
        app.logger.error(f"Error downloading file: {e}")
        return jsonify({'error': 'An error occurred'}), 500

# Route to delete output files
@app.route('/delete_outputs', methods=['DELETE'])
def delete_outputs():
    app.logger.debug('Received request to delete outputs')
    session_id = request.args.get('sessionId')
    if not session_id:
        return jsonify({'error': 'Session ID is required'}), 400

    try:
        session_dir = os.path.join('/tmp', session_id)
        if os.path.exists(session_dir):
            shutil.rmtree(session_dir)
            app.logger.debug(f'Deleted directory: {session_dir}')
        return jsonify({'message': 'Output files deleted successfully'}), 200
    except Exception as e:
        app.logger.error(f"Error deleting output files: {e}")
        return jsonify({'error': 'Error deleting output files'}), 500

##############
### OH YUH ###
##############
# Run the THANG
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8080, debug=True)