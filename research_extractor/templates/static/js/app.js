$(document).ready(function () {

    setupQueryBuilder();
    setupAdditionalFormControls();

    var historyTable = $('#history-table').DataTable({
        responsive: true,
        columns: [
            {
                data: "research_db",
                title: "Datasource"
            },
            {
                data: "search_name",
                title: "Query"
            },
        ],
        select: {
            style: 'single'
        }
    });
    $('#history-table tbody').on('click', 'tr', function () {
        $(this).toggleClass('selected');
    });

    $('#example').DataTable({
        select: {
            style: 'single'
        }
    });

    var resultsTable = $('#results-table').DataTable({
        responsive: true,
        dom: 'Bfrtip',
        buttons: [
            {
                text: '<i class="fa-solid fa-cloud-arrow-up"></i> Curate',
                className: 'buttons-curate',
                action: function (e, dt, button, config) {
                    var data = getTableData();
                    updateSearchResults(data);
                }
            },
            'copy', 'excel', 'pdfHtml5', 'print', 'colvis'
            
        ],
        columns: [
            {
                data: "liked",
                title: 'Like',
                wrap: true,
                render: function (data, type, row) {
                    if (type == 'display') {
                        if (data == true)
                            return '<button class="button button-like liked" value="' + data + '" data-><i class="fa fa-heart"></i><span>Liked</span></button>';
                        else
                            return '<button class="button button-like" value="' + data + '" data-><i class="fa fa-heart"></i><span>Like</span></button>';
                    }
                    return data;
                }
            },
            {
                data: "datasource",
                title: "Data source"
            },
            {
                data: "title",
                title: "Title"
            },
            {
                data: "author",
                title: "author"
            },
            {
                data: "affiliation_country",
                title: "Affiliation Country"
            },
            {
                data: "publication_name",
                title: "Publication Name"
            },
            {
                data: "issn",
                title: "ISSN"
            },
            {
                data: "affiliation_name",
                title: "Affiliation Name"
            },
            {
                data: "url",
                title: "Link",
                render: function (data, type, row) {
                    if (type == 'display') {
                        return '<a target="_blank" href="' + data + '"><i class="fa fa-external-link" aria-hidden="true"></i> Link</a>';
                    }
                    return data;
                }
            }
        ]
    });

    setupDataFramePlugins();

    $('#search').on('click', function () {
        query = $('#query-builder').queryBuilder('getSQL')['sql'];
        query = query.replaceAll('keyword = ', '');
        db = $('#research-db').val();
        force_search = document.getElementById("force-search").checked;
        search_db(db, query, force_search);
    });

    function search_db(db, query, force_search) {
        $('#query-response').html(createQueryText(db, query));
        $.ajax({
            url: "search/research_db/" + db,
            type: 'get',
            beforeSend: () => {
                resultsTable.clear().draw();
                showLoadingScreen(true);
            },
            data: {
                search_text: query,
                force_search: force_search
            },
            success: function (data) {
                resultsTable.rows.add(data.results).draw();
                showLoadingScreen(false);
            },
            error: function (error) {
                alert("Request failed: " + error);
                showLoadingScreen(false);
            }
        });
    }

    $('#lmask').hide();


    function setupQueryBuilder() {
        var options = {
            default_filter: 'keyword',
            filters: [{
                id: 'keyword',
                label: 'Keyword',
                type: 'string',
                size: 200,
                operators: ['equal']
            }]
        };
        $('#query-builder').queryBuilder(options);
    }

    function setupAdditionalFormControls() {
        var search_button_tag = '<button id="search" class="btn btn-md btn-primary pull-right"><span class="glyphicon glyphicon-search"></span> Search</button>'
        var research_db_dropdown_tag = '<label class="form-select" for="research-db">Research DB </label><select id="research-db" class="form-select"></select>'

        $('[data-toggle="tooltip"]').tooltip();   

        $('#query-builder_group_0').append('<div class="search-container">' + research_db_dropdown_tag + search_button_tag + '</div>');

        var research_databases = [
            {
                'db_name': 'scopus',
                'db_desc': 'Scopus/Elsevier'
            },
            {
                'db_name': 'pubmed',
                'db_desc': 'Pubmed',
            },
            {
                'db_name': 'wos',
                'db_desc': 'Web of Science'
            },
            {
                'db_name': 'ieee',
                'db_desc': 'IEEE'
            },
            {
                'db_name': 'all',
                'db_desc': 'All databases'
            }
        ];
        addOptionTags($('#research-db'), research_databases, 'db_name', 'db_desc')
    };

    function addOptionTags(select, options, value_key, label_key) {
        for (const option of options) {
            select.append(createOptionTag(option, value_key, label_key));
        }
    }

    function createOptionTag(option, value_key, label_key) {
        return '<option value="' + option[value_key] + '">' + option[label_key] + '</option>';
    }

    $("#myModal").on("show.bs.modal", function (e) {
        $.ajax({
            url: "search/queries/",
            type: 'get',
            beforeSend: () => {
                historyTable.clear().draw();
                showLoadingScreen(true);
            },
            success: function (results) {
                historyTable.rows.add(results).draw();
                showLoadingScreen(false);
            },
            error: function (error) {
                alert("Request failed: " + error.error);
                console.error(error);
                showLoadingScreen(false);
            }
        });
    });

    $('#search-history').on('click', function () {
        query_details = historyTable.rows('.selected').data()[0];
        if (query_details == undefined)
            alert('Please select a row to load search results!');
        else
            search_db(query_details['research_db'], query_details['search_name'], false);
    });

    $('#results-table tbody').on('click', 'button', function () {
        isLiked = toggleLikeButton($(this));
        resultsTable.cell(this.closest('td')).data(isLiked);
    });

    function toggleLikeButton(button) {
        isLiked = !(button.val() === 'true')
        $(button).val(isLiked);
        return isLiked;
    }

    function setupDataFramePlugins() {
        $('.buttons-curate').removeClass('dt-button').addClass('btn');
        $('.buttons-copy').removeClass('dt-button').addClass('btn btn-primary').html('<i class="fa-regular fa-copy"></i> Copy');
        $('.buttons-excel').removeClass('dt-button').addClass('btn btn-success').html('<i class="fas fa-file-excel"></i> Export');
        $('.buttons-pdf').removeClass('dt-button').addClass('btn btn btn-danger').html('<i class="fa-solid fa-file-pdf"></i> PDF');
        $('.buttons-print').removeClass('dt-button').addClass('btn btn btn-warning').html('<i class="fa-solid fa-print"></i> Print');
    }

    function getTableData() {
        const rows = $('#results-table').DataTable().rows().data();
        var data = [];
        for (var i = 0; i < rows.length; i++) {
            data.push(rows[i]);
        }
        return data;
    }

    function updateSearchResults(data) {
        var queryDetails = getQueryDetails($('#query-response').html());
        request = {
            search_name: queryDetails['search_name'],
            research_db: queryDetails['research_db'],
            results: data
        };
        $.ajax({
            url: "search/update_result/",
            type: 'PUT',
            data: JSON.stringify(request),
            beforeSend: () => showLoadingScreen(true),
            success: function (data) {
                alert('Curation updated!')
                showLoadingScreen(false);
            },
            error: function (error) {
                console.error(error);
                alert('Failed to update curated results!');
                showLoadingScreen(false);
            }
        });
    }

    function getQueryDetails(text) {
        let queryDetails = {};
        let queryContent = text.split(', ');
        queryDetails['research_db'] = queryContent[0].split(': ')[1];
        queryDetails['search_name'] = queryContent[1].split(': ')[1];
        return queryDetails;
    }

    function showLoadingScreen(value) {
        if (value) {
            $('#lmask').show();
        } else {
            $('#lmask').hide();
        }
    }

    function createQueryText(db, query) {
        return "Datasource: " + db + ", Query: " + query;
    }
});


