$(document).ready(function () {

    setupQueryBuilder();
    setupAdditionalFormControls();

    var table = $('#table_id').DataTable({
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
        table.clear().draw();
        // console.log(JSON.stringify($('#query-builder').queryBuilder('getSQL'), undefined, 4));
        query = $('#query-builder').queryBuilder('getSQL')['sql'];
        query = query.replaceAll('keyword = ', '');
        db = $('#research-db').val();
        $('#query-response').html("Datasource: " + db + ", Query: " + query);
        $.ajax({
            url: "search/research_db/" + db,
            type: 'get',
            beforeSend: () => showLoadingScreen(true),
            data: {
                search_text: query,
                force_search: document.getElementById("force-search").checked
            },
            success: function(data) {
                table.rows.add(data.results).draw();
                showLoadingScreen(false);
            },
            error: function(error) {
                alert("Request failed: " + error);
                showLoadingScreen(false);
            }
        });
    });

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
        var $select = $('#research-db')

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
            }
        ];
        for (const research_db of research_databases) {
            tag = '<option value="' + research_db['db_name'] + '">' + research_db['db_desc'] + '</option>';
            $select.append(tag);
        }
    };

    $('#table_id tbody').on('click', 'button', function () {
        isLiked = toggleLikeButton($(this));
        table.cell(this.closest('td')).data(isLiked);
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
        const rows = $('#table_id').DataTable().rows().data();
        var data = [];
        for (var i = 0; i < rows.length; i++) {
            data.push(rows[i]);
        }
        return data;
    }

    function updateSearchResults(data) {
        var queryDetails = getQueryDetails();
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

    function getQueryDetails() {
        let queryDetails = {};
        let queryContent = $('#query-response').html().split(', ');
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
});


