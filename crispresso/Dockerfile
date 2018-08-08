FROM lucapinello/crispresso
# `pip install pysam` installs a broken pysam for some unknown reason :(
RUN conda install -c bioconda pysam
RUN pip install flask boto3
RUN mkdir /code
WORKDIR /code
ADD . /code/

