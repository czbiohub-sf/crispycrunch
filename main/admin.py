from django.contrib import admin

from .models import *
from .forms import *

admin.site.register(Researcher)
admin.site.register(Experiment)
admin.site.register(GuideSelection)
admin.site.register(GuidePlateLayout)
admin.site.register(PrimerDesign)
admin.site.register(PrimerSelection)
admin.site.register(PrimerPlateLayout)


@admin.register(GuideDesign)
class GuideDesignAdmin(admin.ModelAdmin):
    form = GuideDesignForm
