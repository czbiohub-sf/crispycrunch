from abc import abstractmethod
from typing import Tuple

from django.db import models
from .scraperequest import AbstractScrapeRequest


class BatchWebRequest:
    """
    Manages a parallel batch of web requests. Intermediate and final results are
    saved into a passed-in Django model instance.

    There are two public methods: start and get_batch_status.

    There are two abstract methods to override: instance_field and request_key.
    """

    def __init__(
            self,
            web_requester: AbstractScrapeRequest,
            model_instance: models.Model) -> None:
        pass

    @property
    @abstractmethod
    @staticmethod
    def instance_field() -> str:
        """The field name of model_instance that should contain request results."""

    @property
    @abstractmethod
    @staticmethod
    def request_key() -> Tuple[str, ...]:
        """
        A tuple of key names that identify a request result. For example,
        ('target', 'pam_id').
        """
