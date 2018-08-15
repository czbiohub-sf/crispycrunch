from django.urls import path

from main import views


urlpatterns = [
    path('', views.index, name='index'),

    # TODO (gdingle): temp: remove me
    path('crispresso/', views.crispresso, name='crispresso'),

    # BEGIN EXPERIMENT CREATION URLS
    path('experiment/', views.ExperimentView.as_view(), name='Experiment'),

    # Urls here are structured as:
    # The parent object followed by the type of the object to-be-created.
    path('experiment/<int:id>/guide-design/', views.GuideDesignView.as_view(), name='Guide Design'),
    path('guide-design/<int:id>/guide-selection/', views.GuideSelectionView.as_view(), name='Guide Selection'),
    path('guide-design/<int:id>/progress/', views.GuideDesignProgressView.as_view(), name='Guide Design Progress'),
    # TODO (gdingle): do we still want this?
    # path('guide-selection/<int:id>/guide-plate-layout/',
    #      views.GuidePlateLayoutView.as_view(), name='Guide Plate Layout'),

    path('guide-selection/<int:id>/primer-design/', views.PrimerDesignView.as_view(), name='Primer Design'),
    path('primer-design/<int:id>/progress/', views.PrimerDesignProgressView.as_view(), name='Primer Design Progress'),
    path('primer-design/<int:id>/primer-selection/', views.PrimerSelectionView.as_view(), name='Primer Selection'),
    # TODO (gdingle): do we still want this?
    # path('primer-selection/<int:id>/primer-plate-layout/',
    #      views.PrimerPlateLayoutView.as_view(), name='Primer Plate Layout'),
    path('analysis/', views.AnalysisView.as_view(), name='Analysis'),
    # END EXPERIMENT CREATION URLS

    path('primer-selection/<int:id>/experiment-summary/',
         views.ExperimentSummaryView.as_view(), name='Experiment Summary'),
    path('analysis/<int:id>/results/', views.ResultsView.as_view(), name='Experiment Analysis'),

    path('guide-selection/<int:id>/order-form',
         views.GuideOrderFormView.as_view()),
    path('primer-selection/<int:id>/order-form',
         views.PrimerOrderFormView.as_view()),
]
