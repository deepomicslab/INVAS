Bootstrap: docker
From: debian:bullseye-slim

%files
    invas_data_prep.tar.gz /opt/
    invas_assembly.tar.gz /opt/
    scripts/ /app/scripts/
    include/ /app/include/
    src/ /app/src/
    utils/ /app/utils/
    db/ /app/db/
    CMakeLists.txt /app/
    main.cpp /app/
    run_invas.sh /app/

%post
    # 安装编译工具
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        && rm -rf /var/lib/apt/lists/*

    # 解压第一个环境
    mkdir -p /opt/conda/envs/invas_data_prep
    tar -xzf /opt/invas_data_prep.tar.gz -C /opt/conda/envs/invas_data_prep
    rm /opt/invas_data_prep.tar.gz
    /opt/conda/envs/invas_data_prep/bin/conda-unpack

    # 解压第二个环境
    mkdir -p /opt/conda/envs/invas_assembly
    tar -xzf /opt/invas_assembly.tar.gz -C /opt/conda/envs/invas_assembly
    rm /opt/invas_assembly.tar.gz
    /opt/conda/envs/invas_assembly/bin/conda-unpack

    # 克隆并编译 SSW 库
    cd /app
    git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git
    cd Complete-Striped-Smith-Waterman-Library/src
    make
    cp libssw.so /app/scripts/full_pipe/bin/
    cd /app
    rm -rf Complete-Striped-Smith-Waterman-Library

    # 编译 INVAS 核心组件
    export PATH="/opt/conda/envs/invas_assembly/bin:$PATH"
    cd /app/scripts/full_pipe/bin
    cmake ../../../
    make
    chmod +x *

    chmod +x /app/run_invas.sh

%environment
    export PATH="/opt/conda/envs/invas_assembly/bin:$PATH"
    export LC_ALL=C

%runscript
    exec /app/run_invas.sh "$@"

%labels
    Author wangxuedong
    Version 1.0
    Description INVAS - Intragenic Inversion Analysis Software