from django.contrib import admin


from .forms import *
from .models import *


class GuideDesignAdmin(admin.ModelAdmin):

    # The default admin html input errors on ENST00000617316,N,
    # ENST00000278840,C input, so we make read-only.
    # TODO (gdingle): fix with formfield_overrides
    readonly_fields = ['targets_raw']


admin.site.register(Experiment)
admin.site.register(GuideSelection)
admin.site.register(PrimerDesign)
admin.site.register(PrimerSelection)
admin.site.register(Analysis)
admin.site.register(GuideDesign, GuideDesignAdmin)
