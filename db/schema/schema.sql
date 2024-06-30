/* build.sql */
/* Builds the genomicsdb for Luka Culibrk's WGS pipelines */
/* Author: Luka Culibrk */

CREATE DATABASE metadb;

USE metadb;

CREATE TABLE studies(
    id INT unsigned NOT NULL AUTO_INCREMENT,
    rname VARCHAR(150),
    study_pmid VARCHAR(150), 
    ncbi_bioproject_id VARCHAR(50),
    PRIMARY KEY (id),
    UNIQUE KEY unique_rname (rname)
);

CREATE TABLE cells(
    id INT unsigned NOT NULL AUTO_INCREMENT,
    rname VARCHAR(150),
    PRIMARY KEY (id),
    UNIQUE KEY unique_rname (rname)
);

CREATE TABLE patients(
    id INT unsigned NOT NULL AUTO_INCREMENT,
    rname VARCHAR(150),
    germline_sample INT unsigned,
    UNIQUE KEY unique_rname (rname)
);

CREATE TABLE samples(
    id INT unsigned NOT NULL AUTO_INCREMENT,
    rname VARCHAR(150),
    study_id INT unsigned,
    biosample_id VARCHAR(50),
    treatment VARCHAR(1000),
    sample_parent_id INT unsigned,
    cell_id INT unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (study_id)
        REFERENCES studies(id),
    FOREIGN KEY (cell_id)
        REFERENCES cells(id),
    FOREIGN KEY (sample_parent_id)
        REFERENCES samples(id),
    FOREIGN KEY (patients)
        REFERENCES patients(id),
    UNIQUE KEY unique_rname (rname)
);

CREATE TABLE runs(
    id INT unsigned NOT NULL AUTO_INCREMENT,
    rname VARCHAR(150),
    sample_id INT unsigned NOT NULL,
    study_id INT unsigned NOT NULL,
    source VARCHAR(150) NOT NULL,
    ega_id VARCHAR(20),
    sra_id VARCHAR(20),
    local_path VARCHAR(150),
    PRIMARY KEY (id),
    UNIQUE KEY unique_rname (rname),
    FOREIGN KEY (sample_id)
        REFERENCES samples(id),
    FOREIGN KEY (study_id)
        REFERENCES studies(id)
);
CREATE TABLE analyses (
    id INT unsigned NOT NULL AUTO_INCREMENT PRIMARY KEY,
    pipeline_name VARCHAR(20),
    pipeline_version VARCHAR(20),
    analysis_type VARCHAR(20),
    input_table ENUM('studies', 'samples', 'runs') NOT NULL,
    studies_id INT unsigned,
    samples_id INT unsigned,
    runs_id INT unsigned,
    cells_id INT unsigned,
    analysis_time DATETIME DEFAULT CURRENT_TIMESTAMP,
    analysis_dir VARCHAR(250),
    analysis_complete ENUM('False', 'True') DEFAULT 'False',
    FOREIGN KEY (studies_id) REFERENCES studies(id),
    FOREIGN KEY (samples_id) REFERENCES samples(id),
    FOREIGN KEY (runs_id) REFERENCES runs(id),
    FOREIGN KEY (cells_id) REFERENCES cells(id)
);