FROM rhub/r-minimal as base
RUN installr -d -t "curl-dev libxml2-dev libpng-dev gfortran"\
  -a "libcurl libxml2 libpng"\
  DESeq2

RUN installr -d doParallel jsonlite RSQLite parallel

WORKDIR /work

ADD dds2.RDS genesymbol.csv /work 
ADD parallel-raw-data.R merge-db.R /work
CMD ["Rscript", "parallel-raw-data.R"]