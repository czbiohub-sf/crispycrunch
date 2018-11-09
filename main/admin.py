from django.contrib import admin

from .forms import *
from .models import *

# TODO (gdingle): still needed?
# admin.site.register(Researcher)
admin.site.register(Experiment)
admin.site.register(GuideSelection)
admin.site.register(PrimerDesign)
admin.site.register(PrimerSelection)
admin.site.register(Analysis)
admin.site.register(GuideDesign)
