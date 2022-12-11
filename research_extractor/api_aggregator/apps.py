from django.apps import AppConfig


class ApiAggregatorConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "research_extractor.api_aggregator"
    api_key = "b395bcd63c8fb8ee197217a1cb373cef" #35a11e0bb594e5902a369b4cedfbb209
    wos_api_key = "fd0fc81109724d16f0dabcec852249c9232e1d3e"
