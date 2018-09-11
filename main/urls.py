from django.urls import path

from main import views


urlpatterns = [
    path('', views.IndexView.as_view(), name='Crispy Crunch'),

    path('experiment/', views.ExperimentView.as_view(), name='Experiment'),

    # Urls here are structured as:
    # The parent object followed by the type of the object to-be-created.
    path('experiment/<int:id>/guide-design/', views.GuideDesignView.as_view(), name='Guide Design'),
    path('guide-design/<int:id>/progress/', views.GuideDesignProgressView.as_view(), name='Guide Design Progress'),
    path('guide-design/<int:id>/guide-selection/', views.GuideSelectionView.as_view(), name='Guide Selection'),
    path('guide-selection/<int:id>/primer-design/', views.PrimerDesignView.as_view(), name='Primer Design'),
    path('primer-design/<int:id>/progress/', views.PrimerDesignProgressView.as_view(), name='Primer Design Progress'),
    path('primer-design/<int:id>/primer-selection/', views.PrimerSelectionView.as_view(), name='Primer Selection'),

    path('primer-selection/<int:id>/experiment-summary/',
         views.ExperimentSummaryView.as_view(), name='Experiment Summary'),

    path('analysis/', views.AnalysisView.as_view(), name='Analysis'),
    path('analysis/<int:id>/progress/', views.AnalysisProgressView.as_view(), name='Analysis Progress'),
    path('analysis/<int:id>/results/', views.ResultsView.as_view(), name='Results'),

    path('guide-selection/<int:id>/order-form',
         views.GuideOrderFormView.as_view()),
    path('primer-selection/<int:id>/order-form',
         views.PrimerOrderFormView.as_view()),
]
