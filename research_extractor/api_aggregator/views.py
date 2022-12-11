from __future__ import print_function
from django.shortcuts import render
from django.http import JsonResponse
from collections import namedtuple
from rest_framework.parsers import JSONParser
import requests
#scopus
from pybliometrics.scopus.utils import deduplicate, get_freetoread, listify
from django.apps import apps
#pubmed
from pymed import PubMed
from Bio import Entrez
#wos
import time
import woslite_client
from woslite_client.rest import ApiException
from pprint import pprint
import json
from .models import SearchResults #PapersArchive, SearchArchive, FavoritesArchive
from .serializers import Search_ResultsSerializer #SearchArchiveSerializer
import json
import ast

json_parser = JSONParser()

def create_doc():
    fields = 'eid doi pii pubmed_id title subtype subtypeDescription ' \
                'creator afid affilname affiliation_city ' \
                'affiliation_country author_count author_names author_ids url '\
                'author_afids coverDate coverDisplayDate publicationName '\
                'issn source_id eIssn aggregationType volume '\
                'issueIdentifier article_number pageRange description '\
                'authkeywords citedby_count openaccess freetoread '\
                'freetoreadLabel fund_acr fund_no fund_sponsor'
    return namedtuple('Document', fields)


def _join(item, key, sep=";"):
    # """Auxiliary function to join same elements of a list of dictionaries if
    # the elements are not None.
    # """
    try:
        return sep.join([d[key] or "" for d in item["affiliation"]])
    except (KeyError, TypeError):
        return None
    
def _joinurlTo(item, key, sep=";"):
    try:
        for d in item["link"]:
            if d["@ref"]=="scopus":
                return d["@href"]

    except (KeyError, TypeError):
        return None

def _replace_none(lst, repl=""):
    return [repl if v is None else v for v in lst]

def get_app_config():
    return apps.get_app_config('api_aggregator')

def get_api_key():
    return get_app_config().api_key

def search_scopus(request):
    search_text = request.GET.get('search_text')
    print(f"Searching text: {search_text}")
    
    #checking if search alredy present
    searchid = SearchResults.objects.filter(search_name=search_text)
    searchlist = Search_ResultsSerializer(searchid, many=True)

    prevSearch=""
    prevId=""
    prevresults=[]
    # print("searchlist.data     :",searchlist.data)
    if len(searchlist.data)>0:
        for val in searchlist.data:
            # print("seearchname----",val["search_name"])
            prevSearch=val["search_name"]
            prevId=val["search_id"]
            if(len(prevresults)==0):
                prevresults=val["results"]
            else:
                prevresults.extend(val["results"])

    # print("major length------",prevresults)

    #Uncomment start
    api_url = f"https://api.elsevier.com/content/search/scopus?query=all({search_text})&apiKey={get_api_key()}"
    response1 = requests.get(api_url)

    responsearray = response1.json()["search-results"]["entry"]

    docresult = []
    doc = create_doc()
    for item in responsearray:
        info = {}
        # Parse affiliations
        info["affilname"] = _join(item, 'affilname')
        info["afid"] = _join(item, 'afid')
        info["aff_city"] = _join(item, 'affiliation-city')
        info["aff_country"] = _join(item, 'affiliation-country')
        info["url"] = _joinurlTo(item,'link')
        # Parse authors
        try:
            # Deduplicate list of authors
            authors = deduplicate(item['author'])
            # Extract information
            surnames = _replace_none([d['surname'] for d in authors])
            firstnames = _replace_none([d['given-name'] for d in authors])
            info["auth_names"] = ";".join([", ".join([t[0], t[1]]) for t in
                                            zip(surnames, firstnames)])
            info["auth_ids"] = ";".join([d['authid'] for d in authors])
            affs = []
            for auth in authors:
                aff = listify(deduplicate(auth.get('afid', [])))
                affs.append('-'.join([d['$'] for d in aff]))
            if [a for a in affs if a]:
                info["auth_afid"] = ';'.join(affs)
            else:
                info["auth_afid"] = None
        except KeyError:
            pass
        date = item.get('prism:coverDate')
        if isinstance(date, list):
            date = date[0].get('$')
        default = [None, {"$": None}]
        freetoread = get_freetoread(item, ["freetoread", "value"], default)
        freetoreadLabel = get_freetoread(item, ["freetoreadLabel", "value"], default)
        new = doc(article_number=item.get('article-number'),
                    title=item.get('dc:title'), fund_no=item.get('fund-no'),
                    fund_sponsor=item.get('fund-sponsor'),
                    subtype=item.get('subtype'), doi=item.get('prism:doi'),
                    subtypeDescription=item.get('subtypeDescription'),
                    issn=item.get('prism:issn'), creator=item.get('dc:creator'),
                    affilname=info.get("affilname"),
                    author_names=info.get("auth_names"),
                    coverDate=date, volume=item.get('prism:volume'),
                    coverDisplayDate=item.get('prism:coverDisplayDate'),
                    publicationName=item.get('prism:publicationName'),
                    source_id=item.get('source-id'), author_ids=info.get("auth_ids"),
                    aggregationType=item.get('prism:aggregationType'),
                    issueIdentifier=item.get('prism:issueIdentifier'),
                    pageRange=item.get('prism:pageRange'),
                    author_afids=info.get("auth_afid"),
                    affiliation_country=info.get("aff_country"),
                    citedby_count=int(item['citedby-count']),
                    openaccess=int(item['openaccess']),
                    freetoread=freetoread, freetoreadLabel=freetoreadLabel,
                    eIssn=item.get('prism:eIssn'),
                    author_count=item.get('author-count', {}).get('$'),
                    affiliation_city=info.get("aff_city"), afid=info.get("afid"),
                    description=item.get('dc:description'), pii=item.get('pii'),
                    authkeywords=item.get('authkeywords'), eid=item.get('eid'),
                    fund_acr=item.get('fund-acr'), pubmed_id=item.get('pubmed-id'),
                    url=info["url"])
        docresult.append(new)
        
    response = []
    for item in docresult:
        docobj={
                "title": item.title,
                "author":"dummyr",
                "affiliation_country": item.affiliation_country,
                "affiliation_name": item.affilname,
                # "affiliation_country": item.affiliation_country,
                "publication_name": item.publicationName,
                "issn": item.issn,
                "affiliation_name": item.affilname,
                "url": item.url,
                "abstract": '',
                "liked": False
            }

        response.append(
            docobj
        )

    response.extend(prevresults)

    seen = []
    print("prev length",len(response))
    new_list=[]
    for data in response:
        if data["title"] not in seen:
            seen.append(data["title"])
            new_list.append(data)
            
    print("new length",len(new_list))

    # print(len(response))
    searchobj={}
    if prevSearch == "":
        searchobj={
            # "search_id": models.AutoField(primary_key=True)
            "search_name": search_text,
            # search_date = models.DateField(d)
            "research_db": "Scopus",
            "results": new_list
        }
        saveSearch= SearchResults()
        saveSearch.search_name=searchobj['search_name']
        saveSearch.research_db=searchobj["research_db"]
        saveSearch.results=searchobj['results']
        saveSearch.save()

    else:
        searchobj={
            # "search_id": models.AutoField(primary_key=True)
            "search_name": prevSearch,
            # search_date = models.DateField(d)
            "research_db": "Scopus",
            "results": new_list
        }
        # searchid = SearchResults.objects.filter(search_name=search_text)
        # searchlist = Search_ResultsSerializer(searchid, many=True)
        
        oldSR = SearchResults.objects.get(search_id=prevId)
        oldSR.results = searchobj['results']
        oldSR.save() 

    print("returning-----:", searchobj["results"])
    return JsonResponse( searchobj["results"], safe=False)


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
  
def search_pubmed(request):
    search_text = request.GET.get('search_text')
    results = search_query(search_text) #search word = Image Driven 2D Material
    id_list = results['IdList']  # type: ignore
    papers = fetch_details(id_list)
    records=[]
    for i, paper in enumerate(papers['PubmedArticle']):  # type: ignore
        allaffils=''
        name=''
        date=''
        for articledate in paper['MedlineCitation']['Article']['ArticleDate']:
            date=articledate["Year"]+"-"+articledate["Month"]+"-"+articledate["Day"]
        for affillist in paper['MedlineCitation']['Article']['AuthorList']:
            for affil in affillist["AffiliationInfo"]:
                allaffils=affil["Affiliation"]+";"+allaffils
            name=str(affillist["LastName"])+" "+ str(affillist["ForeName"])+"; "+name
        record = {
            'title' : paper['MedlineCitation']['Article']['ArticleTitle'],
            'article_date' : date,
            'creator': name,
            'affiliation_country' : paper['MedlineCitation']['MedlineJournalInfo']['Country'],
            'publicationName' : paper['MedlineCitation']['Article']['Journal']['Title'],
            'issn' : paper['MedlineCitation']['Article']['Journal']['ISSN'],
            'affilname' : allaffils,
            'liked': False
        }
        records.append(record)        
    return JsonResponse(records, safe=False)


# Configure API key authorization: key

def get_api_key():
    return get_app_config().wos_api_key

configuration = woslite_client.Configuration()
configuration.api_key['X-ApiKey'] = get_api_key()

def search_wos(request):
    search_text = request.GET.get('search_text')
    # create an instance of the API class
    # search by specific id
    # integration_api_instance = woslite_client.IntegrationApi(woslite_client.ApiClient(configuration))
    # unique_id = 'WOS:000270372400005'  # str | Primary item(s) id to be searched, ex: WOS:000270372400005. Cannot be null or an empty string. Multiple values are separated by comma.
    search_api_instance = woslite_client.SearchApi(woslite_client.ApiClient(configuration))
    database_id = 'WOS'  # str | Database to search. Must be a valid database ID, one of the following: BCI/BIOABS/BIOSIS/CCC/DCI/DIIDW/MEDLINE/WOK/WOS/ZOOREC. WOK represents all databases.
    usr_query = f'TS='+search_text  # str | User query for requesting data, The query parser will return errors for invalid queries.
    count = 5 # int | Number of records returned in the request
    first_record = 1  # int | Specific record, if any within the result set to return. Cannot be less than 1 and greater than 100000.
    # not required
    lang = 'en'  # str | Language of search. This element can take only one value: en for English. If no language is specified, English is passed by default. (optional)
    sort_field = 'PY+D'  # str | Order by field(s). Field name and order by clause separated by '+', use A for ASC and D for DESC, ex: PY+D. Multiple values are separated by comma. (optional)

    records=[]

    # print(type(records))

    try:
    # Find record(s) by user query
        api_response = search_api_instance.root_get(database_id, usr_query, count, first_record, lang=lang,
                                                                sort_field=sort_field)
        # for more details look at the models
        # firstAuthor = api_response.data[0].author.authors[0]
        # print("Response: ")
        # pprint(type(api_response.data))
        # file = json.dumps(api_response.data, indent=4)
        # with open("api_respone_results.json", "w") as outfile:
        #     outfile.write(file)
        # pprint("First author: " + firstAuthor)
        
        doc_object = {
            'title' : 'Not found',
            'article_date' : 'Not found',
            'author' : 'Not found',
            'affiliation_country' : 'Not found',
            'publicationName' : '',
            'issn' : 'Not found',
            'affilname' : 'Not found'
        }

        for vals in api_response.data:
            if vals.other != '':
                if vals.other.contributor_researcher_id_names != '':
                    if vals.other.contributor_researcher_id_names != None:
                        print("vals--------",vals.other.contributor_researcher_id_names)
                        doc_object['affilname'] = vals.other.contributor_researcher_id_names
            if vals.other != '':
                if vals.other.identifier_issn != '':
                    if vals.other.identifier_issn != None:
                        # print("vals----------",vals.other.identifier_issn)
                        doc_object['issn'] = vals.other.identifier_issn[0]
            if vals.source != '':
                if vals.source.source_title != '':
                    if vals.source.source_title != None:
                        # print("publication name\nvals-------",vals.source.source_title)
                        doc_object['publicationName'] = vals.source.source_title[0]
            if vals.source != '':
                if vals.source.published_biblio_date and vals.source.published_biblio_year != '':
                    if vals.source.published_biblio_date and vals.source.published_biblio_year != None:
                        # print("vals------",vals.source.published_biblio_date,"\t",vals.source.published_biblio_year)
                        doc_object['article_date'] = vals.source.published_biblio_date[0]+" "+vals.source.published_biblio_year[0]
            if vals.author != '':
                if vals.author.authors != '':
                    if vals.author.authors != None:
                        # print("vals------",vals.author.authors)
                        doc_object['author'] = vals.author.authors
            if vals.title != '':
                if vals.title.title != '':
                    if vals.title.title != None:
                        doc_object['title'] = vals.title.title
                        # print("vals----\n",vals.title.title)

            records.append(doc_object)
            # pprint(type(records))

    except ApiException as e:
        print("Exception when calling SearchApi->root_get: %s\\n" % e)


    return JsonResponse(records, safe=False)