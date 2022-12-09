from django.apps import AppConfig


class ScopusConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "research_extractor.scopus"
    api_key = "7f59af901d2d86f78a1fd60c1bf9426a"
