from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.conf import settings


class Command(BaseCommand):

    def handle(self, *args, **options):
        if not User.objects.filter(username="admin").exists():
            admin_email = settings.get('ADMIN_EMAIL', "admin@admin.com")
            User.objects.create_superuser("admin", admin_email, "admin")
