from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),

    # BEGIN EXPERIMENT CREATION URLS
    path('experiment/', views.ExperimentView.as_view(), name='Experiment'),

    # Urls here are structured as:
    # The parent object followed by the type of the object to-be-created.
    path(r'experiment/<int:id>/guide-design/', views.GuideDesignView.as_view(), name='Guide Design'),
    path(r'guide-design/<int:id>/guide-selection/', views.GuideSelectionView.as_view(), name='Guide Selection'),
    path(r'guide-selection/<int:id>/guide-plate-layout/',
         views.GuidePlateLayoutView.as_view(), name='Guide Plate Layout'),

    path(r'guide-selection/<int:id>/primer-design/', views.PrimerDesignView.as_view(), name='Primer Design'),
    path(r'primer-design/<int:id>/primer-selection/', views.PrimerSelectionView.as_view(), name='Primer Selection'),
    path(r'primer-selection/<int:id>/primer-plate-layout/',
         views.PrimerPlateLayoutView.as_view(), name='Primer Plate Layout'),
    # END EXPERIMENT CREATION URLS

    path(r'primer-plate-layout/<int:id>/experiment-summary/',
         views.ExperimentSummaryView.as_view(), name='Experiment Summary'),
    path(r'guide-plate-layout/<int:id>/order-form',
         views.GuideOrderFormView.as_view()),
    path(r'primer-plate-layout/<int:id>/order-form',
         views.PrimerOrderFormView.as_view()),
]
