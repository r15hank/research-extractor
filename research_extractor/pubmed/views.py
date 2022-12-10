from django.shortcuts import render
from django.http import JsonResponse
from rest_framework.parsers import JSONParser
from pymed import PubMed
from Bio import Entrez

json_parser = JSONParser()

# Create your views here.

pubmed = PubMed(tool="MyTool", email="amhatre1@binghamton.edu")

def search_query(query):
    Entrez.email = 'amhatre1@binghamton.edu'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            # retmax='1000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle) 
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'amhatre1@binghamton.edu'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results
  
def search(request):
    search_text = request.GET.get('search_text')
    results = search_query(search_text) #search word = Image Driven 2D Material
    id_list = results['IdList']  # type: ignore
    papers = fetch_details(id_list)
    records=[]
    for i, paper in enumerate(papers['PubmedArticle']):  # type: ignore
        allaffils=''
        name=''
        date=''
        # if (paper['MedlineCitation']['Article']['AuthorList']!= False):
        #     if len(paper['MedlineCitation']['Article']['AuthorList'])>0:
        for articledate in paper['MedlineCitation']['Article']['ArticleDate']:
            date=articledate["Year"]+"-"+articledate["Month"]+"-"+articledate["Day"]
        for affillist in paper['MedlineCitation']['Article']['AuthorList']:
            for affil in affillist["AffiliationInfo"]:
                allaffils=affil["Affiliation"]+";"+allaffils
                # for affilInfo in affil:
                # if affil["Identifier"]=="Affiliation":
                    # allaffils=allaffils+";"+affil["Affiliation"]
            name=str(affillist["LastName"])+" "+ str(affillist["ForeName"])+"; "+name
        # print(paper['MedlineCitation']['Article']['AuthorList'])
        record = {
            'title' : paper['MedlineCitation']['Article']['ArticleTitle'],
            'article_date' : date,
            'creator': name,
            'affiliation_country' : paper['MedlineCitation']['MedlineJournalInfo']['Country'],
            'publicationName' : paper['MedlineCitation']['Article']['Journal']['Title'],
            'issn' : paper['MedlineCitation']['Article']['Journal']['ISSN'],
            'affilname' : allaffils
        }
        
        # temp1 = "{}) {}\n".format(i+1, paper['MedlineCitation']['Article']['ArticleTitle'])
        # print(temp1)
        records.append(record)        

    # with open("pubmed_search_results.json", "w") as outfile:
    #     outfile.write(''.join(temp))
    
    return JsonResponse(records, safe=False)
        
# Pretty print the first paper in full to observe its structure

# json_object=json.dumps(papers['PubmedArticle'][0], indent=2)
# # print(json_object)  # type: ignore
# with open("pubmed_article.json", "w") as outfile:
#     outfile.write(json_object)
#-----------------------------------------------------------------------------------------------
# if __name__ == '__main__':
#     results = search_query('Image Driven 2D Material') #search word = Image Driven 2D Material
#     id_list = results['IdList']  # type: ignore
#     papers = fetch_details(id_list)
#     temp=[]
#     for i, paper in enumerate(papers['PubmedArticle']):  # type: ignore
#         temp1 = "{}) {}\n".format(i+1, paper['MedlineCitation']['Article']['ArticleTitle'])
#         # print(temp1)
#         temp.append(temp1)        
#         # print("{}) {}".format(i+1, paper['MedlineCitation']['Article']['ArticleTitle']))

#     with open("pubmed_search_results.json", "w") as outfile:
#         outfile.write(''.join(temp))