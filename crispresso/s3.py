import boto3
import doctest
import os

from concurrent.futures import ThreadPoolExecutor

PLATE_SIZE = 96
DOWNLOAD_DIR = 'input'
FASTQ_SUFFIX = '.fastq.gz'


def download_fastqs(bucket, prefix, overwrite=False):
    """
    Downloads all fastq files from an s3 folder.

    >>> downloads = download_fastqs('jasonli-bucket', 'JasonHDR/96wp1sorted-fastq/')
    >>> downloads[0].startswith(DOWNLOAD_DIR + '/A1-')
    True
    >>> all(d.endswith(FASTQ_SUFFIX) for d in downloads)
    True
    """
    s3 = boto3.client(
        's3',
        # TODO (gdingle): IMPORTANT! REPLACE AWS CREDENTIALS
        aws_access_key_id='AKIAJAXBOJD6WBX3GXKA',
        aws_secret_access_key='87hRVZo6PeBCIQLtb0EIWOXhgW5vTJbPzGgOWr+n',
    )
    response = s3.list_objects(Bucket=bucket, Prefix=prefix, MaxKeys=1000)
    paths = []
    with ThreadPoolExecutor() as pool:
        for key, size in _get_fastqs(response):
            path = DOWNLOAD_DIR + '/' + key.split('/')[-1]
            paths.append(path)
            if not overwrite and os.path.exists(path) and os.path.getsize(path) == size:
                continue
            else:
                pool.submit(
                    s3.download_file,
                    bucket, key, path,
                )
    return paths


def _get_fastqs(response):
    assert response['IsTruncated'] is False
    fastqs = [(obj['Key'], obj['Size']) for obj in response['Contents']
              if obj['Key'].endswith(FASTQ_SUFFIX)]
    assert len(fastqs) <= PLATE_SIZE * 2, 'Expecting reads of a 96-well plate'
    assert len(fastqs) % 2 == 0, 'Expecting paired reads'
    return fastqs


if __name__ == '__main__':
    # aws s3 ls s3://jasonli-bucket/JasonHDR/96wp1sorted-fastq/
    # download_fastqs('jasonli-bucket', 'JasonHDR/96wp1sorted-fastq/', False)
    doctest.testmod()
