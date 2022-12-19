from rest_framework import serializers
from .models import SearchResults, Results #SearchArchive,FavoritesArchive,PapersArchive

# class SearchArchiveSerializer(serializers.ModelSerializer):
#     class Meta:
#         model =  SearchArchive
#         fields = ('id','searchQuery','searchdbName')

# class FavoritesArchiveSerializer(serializers.ModelSerializer):
#     class Meta:
#         model =  FavoritesArchive
#         fields = ('title','creator','pulblicationName','issn','affilname','urlTo','abstract','searchId')

# class PapersArchive(serializers.ModelSerializer):
#     class Meta:
#         model =  PapersArchive
#         fields = ('title','creator','pulblicationName','issn','affilname','urlTo','abstract','abstract')


class ResultsSerializer(serializers.ModelSerializer):
    class Meta:
        model = Results
        fields = ('title','author','publication_name','article_date','issn','affiliation_name','affiliation_name','affiliation_country','url','abstract','liked')


class Search_ResultsSerializer(serializers.ModelSerializer):
    results = ResultsSerializer(many=True)
    class Meta:
        model =  SearchResults
        fields = ('search_id','search_name','research_db','results')
