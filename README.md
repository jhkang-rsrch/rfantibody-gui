# 🧬 RFantibody GUI

[![GitHub](https://img.shields.io/badge/GitHub-RFantibody-blue?style=flat&logo=github)](https://github.com/RosettaCommons/RFantibody)
[![GitHub](https://img.shields.io/badge/GitHub-AlphaFold3-blue?style=flat&logo=github)](https://github.com/google-deepmind/alphafold3)


A user-friendly graphical interface for RFantibody. This web-based GUI makes it easy to access and utilize the powerful features of RFantibody.

## ✨ Key Features
- 🖥️ Intuitive Web Interface
- 🐳 Easy Installation and Execution with Docker
- 📊 Result Visualization and Analysis Tools

## 🚀 Getting Started

### Prerequisites
- [Docker](https://www.docker.com/get-started)
- [Docker Compose](https://docs.docker.com/compose/install/)

### Installation and Setup

1. Clone the repository
```bash
git clone --recurse-submodules https://github.com/jhkang-rsrch/rfantibody-gui.git
```

2. Navigate to project directory
```bash
cd rfantibody-gui
```

3-1. Prepare RFAntibody
```bash
cd third_party/RFantibody
bash include/download_weights.sh
```

3-2. Prepare AlphaFold3
download DB
```bash
cd ../alphafold3
chmod u+x fetch_databases.sh
./fetch_databases.sh ../af3_data/af3_db
```

3-3. Prepare AlphaFold3 parameter file
From Alphafold3 official [Github](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)

To request access to the AlphaFold 3 model parameters, please complete this [file](https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform). Access will be granted at Google DeepMind’s sole discretion. We(Google Deepmind) will aim to respond to requests within 2–3 business days. You may only use AlphaFold 3 model parameters if received directly from Google. Use is subject to these terms of use.

You will get a zst file: af3.bin.zst

Move this into third_party/af3_data/af3_parameter

4. Launch Docker containers
```bash
cd ../..
# You start here once installed
docker compose up -d
```

5. Waiting for intalling RFAntibody
It may takes over 240s.
```bash
docker logs -f rfantibody-worker
```

6. 🌐 Access
Open your browser and visit:
```
http://localhost:2239
```

7. After use
To turn off,
```bash
docker compose down
```


## 💡 Tips
- Initial startup might take some time due to Docker image downloads
- If you encounter any issues, check the Docker logs:
  ```bash
  docker compose logs
  ```
