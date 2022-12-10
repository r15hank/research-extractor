$(document).ready(function () {

    setupQueryBuilder();
    setupAdditionalFormControls();

    var $table = $('#table_id').DataTable({
        "columns": [
            {
                'data': null,
                title: 'Action',
                wrap: true,
                // render: function () {
                //     return '<div class="btn-group"><span id = heart><i class="fa fa-heart-o" aria-hidden="true" ></i></span></div>'
                // } 
                "render": function (data, type, row) { return '<button class="button button-like"><i class="fa fa-heart"></i><span>Like</span></button>' }
            },
            {
                data: "title",
                title: "Title"
            },
            {
                data: "creator",
                title: "Creator"
            },
            {
                data: "affiliation_country",
                title: "Affiliation Country"
            },
            {
                data: "publicationName",
                title: "Publication Name"
            },
            {
                data: "issn",
                title: "ISSN"
            },
            {
                data: "affilname",
                title: "Affiliation Name"
            }
        ]
    });

    $('#search').on('click', function () {
        // console.log(JSON.stringify($('#query-builder').queryBuilder('getSQL'), undefined, 4));
        query = $('#query-builder').queryBuilder('getSQL')['sql'];
        query = query.replaceAll('keyword = ', '');
        db = $('#research-db').val();
        const url = db + "/search?search_text=" + query;
        $.get(url, function (data, status) {
            $table.rows.add(data).draw();
            setupLikeButton();
            console.log("Data: " + data + "\nStatus: " + status);
        });
    });



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
        $('#query-builder_group_0 .rules-group-header .pull-right').append('<button id="search" class="btn btn-xs btn-primary"><span class="glyphicon glyphicon-search"></span> Search</button>');
        $('#query-builder_group_0 .rules-group-header').append('<select id="research-db" class="form-select"></select>');
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
    }

    function setupLikeButton() {
        $('.button-like').on('click', function () {
            $(this).toggleClass("liked");
        });
    }

});

