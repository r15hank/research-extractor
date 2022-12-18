$(document).ready(function () {

    setupQueryBuilder();
    setupAdditionalFormControls();

    var table = $('#table_id').DataTable({
        responsive: true,
        dom: 'Bfrtip',
        buttons: ['copy', 'excel', 'pdfHtml5', 'print', 'colvis'],
        columns: [
            {
                data: "liked",
                title: 'Like',
                wrap: true,
                render: function (data, type, row) {
                    console.log(data, type);
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
                    console.log(data, type);
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
        $('#query-text').html("Datasource: " + db + ", Query:" + query);
        const url = "search/" + db + "?search_text=" + query;
        $('#lmask').show();
        $.ajax({
            url: url,
            type: 'get',
            data: {
                search_text: query,
                force_search: document.getElementById("force-search").checked
            },
            success: function(data) {
                table.rows.add(data).draw();
                $('#lmask').hide();
            },
            error: function(error) {
                alert("Request failed: " + error);
                $('#lmask').hide();
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
        $('.buttons-copy').removeClass('dt-button').addClass('btn btn-primary').html('<i class="fa-regular fa-copy"></i> Copy');
        $('.buttons-excel').removeClass('dt-button').addClass('btn btn-success').html('<i class="fas fa-file-excel"></i> Export');
        $('.buttons-pdf').removeClass('dt-button').addClass('btn btn btn-danger').html('<i class="fa-solid fa-file-pdf"></i> PDF');
        $('.buttons-print').removeClass('dt-button').addClass('btn btn btn-warning').html('<i class="fa-solid fa-print"></i> Print');
    }
});


