from django.shortcuts import render
from django.http import JsonResponse
from collections import namedtuple
from rest_framework.parsers import JSONParser
import requests
from pybliometrics.scopus.utils import deduplicate, get_freetoread, listify
from django.apps import apps

json_parser = JSONParser()


def create_doc():
    fields = 'eid doi pii pubmed_id title subtype subtypeDescription ' \
                'creator afid affilname affiliation_city ' \
                'affiliation_country author_count author_names author_ids '\
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


def _replace_none(lst, repl=""):
    return [repl if v is None else v for v in lst]

def get_app_config():
    return apps.get_app_config('scopus')

def get_api_key():
    return get_app_config().api_key

def search(request):
    search_text = request.GET.get('search_text')

    api_url = f"https://api.elsevier.com/content/search/scopus?query=all({search_text})&apiKey={get_api_key()}"
    response1 = requests.get(api_url)

    responsearray = response1.json()["search-results"]["entry"];
    docresult = []
    doc = create_doc()
    for item in responsearray:
        info = {}
        # Parse affiliations
        info["affilname"] = _join(item, 'affilname')
        info["afid"] = _join(item, 'afid')
        info["aff_city"] = _join(item, 'affiliation-city')
        info["aff_country"] = _join(item, 'affiliation-country')
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
                    fund_acr=item.get('fund-acr'), pubmed_id=item.get('pubmed-id'))
        docresult.append(new)
        response = []
        for item in docresult:
            response.append(
                {
                    "title": item.title,
                    "creator": item.creator,
                    "affiliation_country": item.affiliation_country,
                    "publicationName": item.publicationName,
                    "issn": item.issn,
                    "affilname": item.affilname
                }
            )
    return JsonResponse(response, safe=False)

