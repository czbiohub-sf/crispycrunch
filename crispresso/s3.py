import boto3  # type: ignore # noqa
import doctest
import os

from botocore.config import Config # type: ignore

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List

PLATE_SIZE = 96
BASE_DIR = os.path.dirname(__file__)
DOWNLOAD_DIR = os.path.join(BASE_DIR,
                            # Use git-committed dir in prod to avoid permission issues
                            'fastqs' if 'RDS_DB_NAME' in os.environ else 'input'
                            )
FASTQ_SUFFIX = '.fastq.gz'


def download_fastqs(bucket: str, prefix: str, overwrite=False) -> List[str]:
    """
    Downloads all fastq files from an s3 folder.
    >>> downloads = download_fastqs('jasonli-bucket', 'JasonHDR/96wp1sorted-fastq/')
    >>> downloads[0].startswith(DOWNLOAD_DIR)
    True
    >>> all(d.endswith(FASTQ_SUFFIX) for d in downloads)
    True
    """

    s3 = boto3.client(
        's3',
        config=Config(**{
            'connect_timeout': 5,
            'read_timeout': 300}),
        # See https://boto3.amazonaws.com/v1/documentation/api/latest/guide/configuration.html
        # aws_access_key_id='',
        # aws_secret_access_key='',
    )
    response = s3.list_objects(Bucket=bucket, Prefix=prefix, MaxKeys=1000)
    paths = []
    # Default size of boto thread pool is 10
    # Avoid "Connection pool is full, discarding connection"
    with ThreadPoolExecutor(10) as pool:
        for key, size in _get_fastqs(response):
            new_dir = Path(DOWNLOAD_DIR) / Path(key).parent
            new_dir.mkdir(parents=True, exist_ok=True)
            new_filepath = new_dir / Path(key).name
            paths.append(new_filepath)
            if not overwrite and new_filepath.exists() and new_filepath.stat().st_size == size:
                continue
            else:
                pool.submit(
                    # TODO (gdingle): create new client? see
                    # https://boto3.amazonaws.com/v1/documentation/api/latest/guide/resources.html#multithreading
                    s3.download_file,
                    bucket, key, str(new_filepath),
                )
    # Callers expect strings
    return [str(p) for p in paths]


def _get_fastqs(response) -> list:
    assert response['IsTruncated'] is False

    if 'Contents' not in response:
        raise ValueError('No contents in AWS bucket prefix {}'.format(response['Prefix']))

    fastqs = [(obj['Key'], obj['Size']) for obj in response['Contents']
              if obj['Key'].endswith(FASTQ_SUFFIX)]
    # TODO (gdingle): what should the min and max in a dir be?
    # assert len(fastqs) <= PLATE_SIZE * 2, 'Expecting reads of a 96-well plate'
    # assert len(fastqs) % 2 == 0, 'Expecting paired reads'
    return fastqs


if __name__ == '__main__':
    # aws s3 ls s3://jasonli-bucket/JasonHDR/96wp1sorted-fastq/
    # print(download_fastqs('jasonli-bucket', 'JasonHDR/96wp1sorted-fastq/', False))
    # print(download_fastqs('jasonli-bucket', 'CrispyCrunch/', False))
    # print(download_fastqs('donor6-1', 'fastq_fildes/', False))
    doctest.testmod()
