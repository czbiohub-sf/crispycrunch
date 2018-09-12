"""crispycrunch URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
"""
from django.contrib import admin
from django.urls import include, path
from django.views.generic import RedirectView

urlpatterns = [
    path('admin/', admin.site.urls),
    path('main/', include('main.urls')),

    # See https://stackoverflow.com/questions/9371378/warning-not-found-favicon-ico
    path('favicon.ico', RedirectView.as_view(url='/static/crispy-crunch-logo-small.png')),
]
