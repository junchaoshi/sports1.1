# SPORTS1.1 Demo Data Analysis Guide

## 1. Download Mouse Sperm Small RNA-Seq Demo Data

Use the following command to download the demo data to your desired path:

```bash
fastq-dump --outdir /path_to_your_test_folder/fq_data/ SRR1538547
```

> The download may take approximately 5 minutes, depending on network speed.

---

## 2. Install SPORTS1.1 and Its Dependencies

Follow the steps below to install **SPORTS1.1** and all necessary software.

### Download SPORTS1.1 Pipeline Package

```bash
wget -P /path_to_software_folder/ https://github.com/junchaoshi/SPORTS1.1/archive/master.zip
```

### Download Necessary Software, Packages, and Reference Databases

1. **Perl 5** (Tested versions: v5.14.2, v5.22.1)  
   - [https://www.perl.org](https://www.perl.org)  
   > *Note*: Perl 5 might already be installed on your system.

2. **Bowtie** (Tested versions: 1.1.2, 1.2.1.1)  
   - [http://bowtie-bio.sourceforge.net/index.shtml](http://bowtie-bio.sourceforge.net/index.shtml)

3. **SRA Toolkit** (Tested version: 2.8.2)  
   - [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)

4. **cutadapt** (Tested versions: 1.11, 5.0)  
   - [http://cutadapt.readthedocs.io/en/stable/index.html](http://cutadapt.readthedocs.io/en/stable/index.html)

5. **R** (Tested versions: 3.2.3, 3.2.5)  
   - [https://www.r-project.org/](https://www.r-project.org/)

6. **Reference Database**  
   - See pre-compiled databases on the SPORTS1.1 GitHub page: [https://github.com/junchaoshi/SPORTS1.1/](https://github.com/junchaoshi/SPORTS1.1/)

---

### Installation Steps for SPORTS1.1 and Dependencies

#### Install SPORTS1.1

1. **Unpack SPORTS1.1 package:**

   ```bash
   unzip /path_to_SPORTS1.1_folder/SPORTS1.1-master.zip
   ```

2. **Add SPORTS1.1 to your `PATH`:**

   ```bash
   echo 'export PATH=$PATH:/path_to_software_folder/sports1.1-master/source' >> ~/.bashrc
   chmod 755 /path_to_software_folder/sports1.1-master/source/sports.pl
   ```

#### Install Bowtie

1. **Unpack Bowtie:**

   ```bash
   unzip /path_to_software_folder/bowtie-1.x.x-linux-x86_64.zip
   ```

2. **Add Bowtie to your `PATH`:**

   ```bash
   echo 'export PATH=$PATH:/path_to_bowtie' >> ~/.bashrc
   ```

3. **(Optional) Install Bowtie as administrator:**

   ```bash
   sudo apt-get install bowtie
   ```

#### Install SRA Toolkit

1. **Unpack SRA Toolkit files.**

2. **Add SRA Toolkit to your `PATH`:**

   ```bash
   echo 'export PATH=$PATH:/path_to_sra-toolkit/bin' >> ~/.bashrc
   ```

#### Install cutadapt

1. **Install the latest version of cutadapt:**

   ```bash
   pip install --user --upgrade cutadapt
   ```

2. **Add cutadapt to your `PATH`:**

   ```bash
   echo 'export PATH=$PATH:$HOME/.local/bin' >> ~/.bashrc
   ```

#### Install R and Required Packages

1. **Unpack the R source file:**

   ```bash
   tar -xf R-x.y.z.tar.gz
   ```

2. **Navigate to the R directory:**

   ```bash
   cd R-x.y.z
   ```

3. **Build and install R:**

   ```bash
   ./configure
   make
   make check
   make install
   ```

4. **Install R packages:**

   ```r
   R
   install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')
   install.packages('data.table', dependencies=TRUE, repos='http://cran.rstudio.com/')
   install.packages('stringr', dependencies=TRUE, repos='http://cran.rstudio.com/')
   q()
   n
   ```

#### Apply Environment Changes

```bash
source ~/.bashrc
```

#### Test Installation

Verify installations by running:

```bash
perl -v
sports.pl -h
bowtie
fastq-dump
cutadapt -h
R --version
```

> If errors occur, reinstall the respective software.

---

## 3. Download Mouse Reference Database

Use the following commands to download and prepare the reference database:

```bash
mkdir -p /path_to_your_test_folder/reference/
wget -P /path_to_your_test_folder/reference/ --no-check-certificate "https://ncrnainfo-my.sharepoint.com/personal/sports_ncrna_info/_layouts/15/guestaccess.aspx?docid=0833653a140eb47f098267d7a23d3b63c&authkey=Ab8aoYC8paqFI2yRabIo7Ok"
mv /path_to_your_test_folder/reference/'guestaccess.aspx?docid=0833653a140eb47f098267d7a23d3b63c&authkey=Ab8aoYC8paqFI2yRabIo7Ok' /path_to_your_test_folder/reference/Mus_musculus.tar.gz
tar -zxvf /path_to_your_test_folder/reference/Mus_musculus.tar.gz
```

---

## 4. Run SPORTS1.1

Execute the following script:

```bash
nohup sports.pl \
   -i /path_to_your_test_folder/fq_data/SRR1538547.fastq \
   -p 4 \
   -M 1 \
   -a -y TGGAATTCTCGGGTGCCAAGGAACTCCAGT \
   -g /path_to_your_test_folder/reference/Mus_musculus/genome/mm10/genome \
   -m /path_to_your_test_folder/reference/Mus_musculus/miRBase/21/miRBase_21-mmu \
   -r /path_to_your_test_folder/reference/Mus_musculus/rRNAdb/mouse_rRNA \
   -t /path_to_your_test_folder/reference/Mus_musculus/GtRNAdb/mm10/mm10-tRNAs \
   -e /path_to_your_test_folder/reference/Mus_musculus/Ensembl/release-89/Mus_musculus.GRCm38.ncrna \
   -f /path_to_your_test_folder/reference/Mus_musculus/Rfam/12.3/Rfam-12.3-mouse \
   -w /path_to_your_test_folder/reference/Mus_musculus/piRBase/piR_mouse \
   -o /path_to_your_test_folder/Demo_SPORTS_output/ &
```

> Running the SPORTS1.1 script for the first time may take about 40 minutes to generate the Bowtie index using 4 threads. Subsequent runs for the sample will take approximately 3 minutes.
