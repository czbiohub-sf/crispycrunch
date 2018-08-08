FROM lucapinello/crispresso
RUN pip install flask boto3
# `pip install pysam` installs a broken pysam for some unknown reason :(
RUN conda install -c bioconda pysam
RUN mkdir /code
WORKDIR /code
ADD . /code/

