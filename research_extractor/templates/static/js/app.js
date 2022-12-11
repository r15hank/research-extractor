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
                title: "Link"
            }
        ]
    });

    $('#search').on('click', function () {
        table.clear().draw();
        // console.log(JSON.stringify($('#query-builder').queryBuilder('getSQL'), undefined, 4));
        query = $('#query-builder').queryBuilder('getSQL')['sql'];
        query = query.replaceAll('keyword = ', '');
        db = $('#research-db').val();
        const url = "search/" + db + "?search_text=" + query;
        $('#lmask').show();
        $.get(url, function (data, status) {
            table.rows.add(data).draw();
            $('#lmask').hide();
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
});


