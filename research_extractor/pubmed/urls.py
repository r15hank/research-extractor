from django.urls import path
from research_extractor.pubmed import views

urlpatterns = [
    path("search/", views.search)
]