from django.core.mail import send_mail
from django.core.management.base import BaseCommand
from django.conf import settings


class Command(BaseCommand):

    def handle(self, *args, **options):
        # See https://docs.djangoproject.com/en/2.1/topics/email/
        send_mail(
            'Subject here',
            'Here is the message.',
            settings.SERVER_EMAIL,
            [settings.ADMIN_EMAIL],
            fail_silently=False,
        )
