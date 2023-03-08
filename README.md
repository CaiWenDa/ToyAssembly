# alphaHiASM

    正在开发中的基因组装工具

## 依赖项

该项目依赖于 boost 和 seqan 库开发，安装前需要先下载它们。

```
seqan: http://packages.seqan.de/seqan-library/seqan-library-2.4.0.zip
boost: https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz
```

## 安装

克隆存储库

```
git clone https://github.com/CaiWenDa/alphaHiASM.git
```

下载 seqan 和 boost 库

```
wget http://packages.seqan.de/seqan-library/seqan-library-2.4.0.zip
wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz
```

创建 seqan 目录

```
mkdir alphaHiASM/seqan
```

解压缩

```
unzip seqan-library-2.4.0.zip
```

将 “seqan-library-2.4.0/*” 复制到 alphaHiASM/seqan

```
cp -r seqan-library-2.4.0/* alphaHiASM/seqan
```

创建 boost 目录

```
mkdir alphaHiASM/boost
```

解压缩

```
tar -zxvf boost_1_78_0.tar.gz
```

将“boost_1_78_0/*”复制到alphaHiASM/boost

```
cp -r boost_1_78_0/* alphaHiASM/boost
```

切换到 alphaHiASM 目录

```
cd alphaHiASM
```

创建 build 目录

```
mkdir build
```

切换到 build 目录

```
cd build
```

编译

```
cmake ../
make
```

尝试运行

```
./alphaHiASM
```

## 使用方法

    alphaHiASM
        -f --file inFileName
    输入文件，格式为 fasta 或 fastq
        -o --output outFileName
    输出组装文件地址，格式为文本文件
        --genomeSize SIZE
    要组装物种的基因组大小，单位为碱基数
        -t --threads int
    线程数，默认为 1
        [--minOverlap SIZE]
    最小重叠长度，默认 2000
        [--help]
    获得帮助信息
        [--paf inOverlapFile]
    读取已经存在的比对文件信息
        [--overlap outOverlapFile]
    输出的比对信息文件名
        [--kmerLen kmer len]
    k-mer 长度，默认为 31
        [--step detect step]
    检测步长，默认为 1

## 示例

将原始读数文件 reads.fasta 进行自比对，比对结果输出到 output.csv 中：

```
alphaHiASM -f reads.fasta --overlap output.csv
```

使用多线程提高效率，以 12 线程为例：

```
alphaHiASM -f reads.fasta --overlap output.csv -t 12
```

将原始读数文件 reads.fasta 进行自比对同时进行组装，组装物种的基因组大小 genomeSize 估计值为 4.6m，比对和组装结果将分别输出到 output.csv 和 asm.fasta 文件中

```
alphaHiASM -f reads.fasta --overlap output.csv --genomeSize 4.6m -o asm.fasta
```

可同样使用多线程提高效率，以 12 线程为例：

```
alphaHiASM -f reads.fasta --overlap output.csv --genomeSize 4.6m -o asm.fasta -t 12
```

## 输出格式

alphaHiASM 的重叠检测结果以文本文件输出，每一行为一条检测距离，格式如表所示
|列 |类型 |描述|
|:---:|:---:|:---:|
|1 |int |查询序列编号|
|2 |int |目标序列编号|
|3 |bool |0 正向比对，1 反向比对|
|4 |int |查询序列开始位置（从0开始）|
|5 |int |查询序列结束位置（从0开始）|
|6 |int |目标序列开始位置（从0开始）|
|7 |int |目标序列结束位置（从0开始）|
|8 |int |查询序列长度|
|9 |int |目标序列长度|
