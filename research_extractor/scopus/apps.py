from django.apps import AppConfig


class ScopusConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "research_extractor.scopus"
    api_key = "b395bcd63c8fb8ee197217a1cb373cef"
