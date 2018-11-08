"""crispycrunch URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
"""
from django.contrib import admin
from django.urls import include, path
from django.views.generic import RedirectView, TemplateView

urlpatterns = [
    path('', TemplateView.as_view(template_name='home.html'), name='CrispyCrunch'),
    path('admin/', admin.site.urls),
    path('main/', include('main.urls')),

    # See https://stackoverflow.com/questions/9371378/warning-not-found-favicon-ico
    path('favicon.ico', RedirectView.as_view(url='/static/biohub-icon.png')),
]

# TODO (gdingle): 500 error pass in exception message with custom view see https://docs.djangoproject.com/en/2.0/ref/views/
