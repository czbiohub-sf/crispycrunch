from django.contrib.auth.decorators import login_required
from django.urls import path

from main import views


urlpatterns = [
    path('',
         login_required(views.IndexView.as_view()),
         name='Index'
         ),

    path('experiment/',
         login_required(views.ExperimentView.as_view()),
         name='Experiment'
         ),
    path('experiment/<int:id>/delete',
         login_required(views.ExperimentDeleteView.as_view())),

    path('signup/', views.CreateUserView.as_view(), name='Sign-up'),

    # Urls here are structured as:
    # The parent object followed by the type of the object to-be-created.
    path('experiment/<int:id>/guide-design/',
         login_required(views.GuideDesignView.as_view()),
         name='Guide Design'
         ),
    path('guide-design/<int:id>/progress/',
         login_required(views.GuideDesignProgressView.as_view()),
         name='Guide Design Progress'
         ),
    path('guide-design/<int:id>/guide-selection/',
         login_required(views.GuideSelectionView.as_view()),
         name='Guide Selection'
         ),
    path('guide-selection/<int:id>/primer-design/',
         login_required(views.PrimerDesignView.as_view()),
         name='Primer Design'
         ),
    path('primer-design/<int:id>/progress/',
         login_required(views.PrimerDesignProgressView.as_view()),
         name='Primer Design Progress'
         ),
    path('primer-design/<int:id>/primer-selection/',
         login_required(views.PrimerSelectionView.as_view()),
         name='Primer Selection'
         ),

    path('primer-selection/<int:primer_selection_id>/experiment-summary/',
         login_required(views.ExperimentSummaryView.as_view()),
         name='Experiment Summary'
         ),
    # Alternate path to summary
    path('experiment/<int:experiment_id>/summary/',
         login_required(
             views.ExperimentSummaryView.as_view()),
         name='Experiment Summary'
         ),

    path(
        'analysis/',
        login_required(views.AnalysisView.as_view()),
        name='Analysis'
    ),
    path(
        'analysis/<int:id>/progress/',
        login_required(views.AnalysisProgressView.as_view()),
        name='Analysis Progress'
    ),
    path(
        'analysis/<int:id>/results/',
        login_required(views.ResultsView.as_view()),
        name='Results'
    ),
    path(
        'analysis/<int:id>/custom/',
        login_required(views.CustomAnalysisView.as_view()),
        name='Custom Analysis'
    ),
    path('analysis/<int:id>/delete',
         login_required(views.AnalysisDeleteView.as_view())),

    # Downloads
    path('guide-selection/<int:id>/order-form',
         login_required(views.GuideOrderFormView.as_view())),
    path('primer-selection/<int:id>/order-form',
         login_required(views.PrimerOrderFormView.as_view())),
    path('primer-selection/<int:id>/illumina-sheet',
         login_required(views.IlluminaSheetView.as_view())),
    path('primer-selection/<int:id>/hdr-order-form',
         login_required(views.UltramerOrderFormView.as_view())),
]
