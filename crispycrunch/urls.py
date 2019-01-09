"""crispycrunch URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
"""
from django.conf import settings
from django.contrib import admin
from django.urls import include, path
from django.views.generic import RedirectView, TemplateView

urlpatterns = [
    path('', TemplateView.as_view(template_name='home.html'), name='CrispyCrunch'),
    path('howto/', TemplateView.as_view(template_name='howto.html'), name='How to use CrispyCrunch'),
    # See LOGIN_REDIRECT_URL
    path('accounts/', include('django.contrib.auth.urls')),
    path('admin/', admin.site.urls),
    path('main/', include('main.urls')),
    # See https://stackoverflow.com/questions/9371378/warning-not-found-favicon-ico
    path('favicon.ico', RedirectView.as_view(url='/static/biohub-icon.png')),
]

if settings.DEBUG:
    import debug_toolbar
    urlpatterns += [
        # See https://django-debug-toolbar.readthedocs.io/en/latest/installation.html
        path('__debug__/', include(debug_toolbar.urls)),
    ]

# TODO (gdingle): 500 error pass in exception message with custom view see https://docs.djangoproject.com/en/2.0/ref/views/
