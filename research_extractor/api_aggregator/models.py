from djongo import models
# Create your models here.
# class PapersArchive(models.Model):
#     title=models.CharField(primary_key=True,max_length=1000)
#     creator=models.CharField(max_length=1000)
#     pulblicationName=models.CharField(max_length=1000)
#     issn=models.CharField(max_length=1000)
#     affilname=models.CharField(max_length=1000)
#     urlTo=models.CharField(max_length=1000)
#     abstract=models.CharField(max_length=1000)
#     abstract=models.CharField(max_length=1000)

# class SearchArchive(models.Model):
#     # searchId=models.IntegerField(primary_key=True)
#     id=models.AutoField(primary_key=True)   
#     searchQuery=models.CharField(max_length=1000)
#     searchdbName=models.CharField(max_length=1000)

# class FavoritesArchive(models.Model):
#     title=models.CharField(primary_key=True,max_length=1000)
#     creator=models.CharField(max_length=1000)
#     pulblicationName=models.CharField(max_length=1000)
#     issn=models.CharField(max_length=1000)
#     affilname=models.CharField(max_length=1000)
#     urlTo=models.CharField(max_length=1000)
#     abstract=models.CharField(max_length=1000)
#     searchId=models.CharField(max_length=1000)


class Results(models.Model):
    title = models.CharField(primary_key=True, max_length=1000)
    author = models.CharField(max_length=1000)
    publication_name = models.CharField(max_length=1000, blank=True, null=True)
    issn = models.CharField(max_length=1000, blank=True, null=True)
    affiliation_name = models.CharField(max_length=2000, blank=True, null=True)
    affiliation_country = models.CharField(max_length=100, blank=True, null=True)
    abstract = models.CharField(max_length=10000, blank=True, null=True)
    url = models.CharField(max_length=1000)
    liked = models.BooleanField(default=False)
    class Meta:
        managed = False


class SearchResults(models.Model):
    search_id = models.AutoField(primary_key=True)
    search_name = models.CharField(max_length=10000)
    # search_date = models.DateField(d)
    research_db = models.CharField(max_length=50)
    # results = models.ArrayReferenceField(to=Results,on_delete=models.CASCADE)#,model_container=Results

    results = models.ArrayField(model_container=Results)
