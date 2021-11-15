/**
 * (c) 2013 Jexcel Plugin v2.1.0 | Bossanova UI
 * http://www.github.com/paulhodel/jexcel
 *
 * @author: Paul Hodel <paul.hodel@gmail.com>
 * @description: Create light embedded spreadsheets on your webpages
 * 
 * ROADMAP:
 * Online collaboration
 */

(function( $ ){

var methods = {

    /**
     * Innitialization, configuration and loading
     * 
     * @param {Object} options configuration
     * @return void
     */
    init : function( options ) {
        // Loading default configuration
        var defaults = {
            // Column types and configurations
            columns:[],
            // Column header titles
            colHeaders:[],
            // Column width sizes
            colWidths:[],
            // Column alignment
            colAlignments:[],
            // Colum header classes
            colHeaderClasses:[],
            // Column width that is used by default
            defaultColWidth:50,
            // Minimal number of blank rows in the end
            minSpareRows:0,
            // Minimal number of blank cols in the end
            minSpareCols:0,
            // Minimal table dimensions
            minDimensions:[0,0],
            // Custom context menu
            contextMenu:null,
            // Allow column sorting
            columnSorting:true,
            // Allow column resizing
            columnResize:true,
            // Allow row dragging
            rowDrag:true,
            // Allow table edition
            editable:true,
            // Allow new rows
            allowInsertRow:true,
            // Allow new rows
            allowManualInsertRow:true,
            // Allow new columns
            allowInsertColumn:true,
            // Allow new rows
            allowManualInsertColumn:true,
            // Allow row delete
            allowDeleteRow:true,
            // Allow column delete
            allowDeleteColumn:true,
            // Allow cell commments
            allowComments:true,
            // Global wrap
            wordWrap:false,
            // Filename
            csvFileName:'jexcel',
            // Disable corner selection
            selectionCopy:true,
            // Allow Overflow
            tableOverflow:false,
            // Allow Overflow
            tableHeight:200,
            // Toolbar
            toolbar:null,
            // About message
            about:'jExcel Spreadsheet\\nVersion 2.1.0\\nAuthor: Paul Hodel <paul.hodel@gmail.com>\\nWebsite: https://bossanova.uk/jexcel'
        };

        // Id
        var id = $(this).prop('id');

        // Configuration holder
        var options =  $.extend(defaults, options);

        // Save the div that will hold the grid
        options.div = this;

        // Compatibility
        if (options.manualColumnResize != undefined) {
            options.columnResize = options.manualColumnResize;
        }
        if (options.manualRowMove != undefined) {
            options.rowDrag = options.manualRowMove;
        }

        // Register options
        if (! $.fn.jexcel.defaults) {
            // Events control
            $.fn.jexcel.ignoreEvents = false;
            $.fn.jexcel.ignoreHistory = false;

            // Configuration container
            $.fn.jexcel.defaults = [];

            // Formula element factory
            $.fn.jexcel.factory = function(k, v) {
                this.k = k;
                this.v = v;
            }

            $.fn.jexcel.factory.prototype.toString = function() {
                if ($.fn.jexcel.current) {
                    return $.fn.jexcel.defaults[$.fn.jexcel.current].values[this.k];
                } else {
                    return this.v;
                }
            }

            // Copy helpers
            $.fn.jexcel.hash = null;
            $.fn.jexcel.style = null;
        }

        // Global configuration
        $.fn.jexcel.defaults[id] = options;

        // Data holder cannot be blank
        if (! $.fn.jexcel.defaults[id].data) {
            $.fn.jexcel.defaults[id].data = [];
        }

        // Min data length
        if (! $.fn.jexcel.defaults[id].data.length) {
            $.fn.jexcel.defaults[id].data = [[]];
        }

        // Create history track array
        $.fn.jexcel.defaults[id].history = [];
        $.fn.jexcel.defaults[id].historyIndex = -1;


        // Load the table data based on an CSV file
        if ($.fn.jexcel.defaults[id].csv) {
            if (! $.csv) {
                // Required lib not present
                console.error('Jexcel error: jquery-csv library not loaded');
            } else {
                // Comma as default
                $.fn.jexcel.defaults[id].delimiter = $.fn.jexcel.defaults[id].delimiter || ',';

                // Load CSV file
                $.ajax({
                    id: id,
                    url: $.fn.jexcel.defaults[id].csv,
                    success: function (result) {
                        // Convert data
                        var data = $.csv.toArrays(result);

                        // Headers
                        if ($.fn.jexcel.defaults[this.id].csvHeaders == true) {
                            $.fn.jexcel.defaults[this.id].colHeaders = data.shift();
                        }

                        // Data
                        $.fn.jexcel.defaults[this.id].data = data;
 
                        // Prepare table
                        $.fn.jexcel('prepareTable', this.id);
                    }
                });
            }
        } else if (options.url) {
            // Load json external file
            $.ajax({
                id: id,
                url: $.fn.jexcel.defaults[id].url,
                dataType:'json',
                success: function (result) {
                    // Data
                    $.fn.jexcel.defaults[this.id].data = (result.data) ? result.data : result;
                     // Prepare table
                    $.fn.jexcel('prepareTable', this.id);
                }
            });
        } else {
            // Prepare table
            $.fn.jexcel('prepareTable', id);
        }
    },

    /**
     * Prepare table
     */
    prepareTable : function (id) {
        // Loading initial data from remote sources
        var results = [];

        // Number of columns
        var size = $.fn.jexcel.defaults[id].colHeaders.length;

        if ($.fn.jexcel.defaults[id].data[0].length > size) {
            size = $.fn.jexcel.defaults[id].data[0].length;
        }

        // Minimal dimensions
        if ($.fn.jexcel.defaults[id].minDimensions[0] > size) {
            size = $.fn.jexcel.defaults[id].minDimensions[0];
        }

        // Preparations
        for (var i = 0; i < size; i++) {
            // Default headers
            if (! $.fn.jexcel.defaults[id].colHeaders[i]) {
                $.fn.jexcel.defaults[id].colHeaders[i] = '';
            }
            // Default column description
            if (! $.fn.jexcel.defaults[id].columns[i]) {
                $.fn.jexcel.defaults[id].columns[i] = { type:'text' };
            } else if (! $.fn.jexcel.defaults[id].columns[i]) {
                $.fn.jexcel.defaults[id].columns[i].type = 'text';
            }
            if (! $.fn.jexcel.defaults[id].columns[i].source) {
                $.fn.jexcel.defaults[id].columns[i].source = [];
            }
            if (! $.fn.jexcel.defaults[id].columns[i].options) {
                $.fn.jexcel.defaults[id].columns[i].options = [];
            }
            if (! $.fn.jexcel.defaults[id].columns[i].editor) {
                $.fn.jexcel.defaults[id].columns[i].editor = null;
            }
            if (! $.fn.jexcel.defaults[id].columns[i].allowEmpty) {
                $.fn.jexcel.defaults[id].columns[i].allowEmpty = false;
            }
            if (! $.fn.jexcel.defaults[id].colWidths[i]) {
                $.fn.jexcel.defaults[id].colWidths[i] = $.fn.jexcel.defaults[id].defaultColWidth || '50';
            }
            if (! $.fn.jexcel.defaults[id].colAlignments[i]) {
                $.fn.jexcel.defaults[id].colAlignments[i] = 'center';
            }
            if (! $.fn.jexcel.defaults[id].colHeaderClasses[i]) {
                $.fn.jexcel.defaults[id].colHeaderClasses[i] = '';
            }

            // Pre-load initial source for json autocomplete
            if ($.fn.jexcel.defaults[id].columns[i].type == 'autocomplete' || $.fn.jexcel.defaults[id].columns[i].type == 'dropdown') {
                // if remote content
                if ($.fn.jexcel.defaults[id].columns[i].url) {
                    results.push($.ajax({
                        id:id,
                        url: $.fn.jexcel.defaults[id].columns[i].url,
                        index: i,
                        dataType:'json',
                        success: function (result) {
                            // Create the dynamic sources
                            $.fn.jexcel.defaults[this.id].columns[this.index].source = result;
                        },
                        error: function (result) {
                            console.error('It was not possible to load the url: ' + this.url);
                        }
                    }));
                }
            } else if ($.fn.jexcel.defaults[id].columns[i].type == 'calendar') {
                // Default format for date columns
                if (! $.fn.jexcel.defaults[id].columns[i].options.format) {
                    $.fn.jexcel.defaults[id].columns[i].options.format = 'DD/MM/YYYY';
                }
            }
        }

        // In case there are external json to be loaded before create the table
        if (results.length > 0) {
            // Waiting all external data is loaded
            $.when.apply(this, results).done(function() {
                // Create the table
                $($.fn.jexcel.defaults[id].div).jexcel('createTable');
            });
        } else {
            // No external data to be loaded, just created the table
            $($.fn.jexcel.defaults[id].div).jexcel('createTable');
        }
    },

    /**
     * Create the table
     * 
     * @return void
     */
    createTable : function() {
        // Id
        var id = $(this).prop('id');

        // Var options
        var options = $.fn.jexcel.defaults[id];

        // Element
        $(this).prop('class', 'jexcel');

        // Create toolbar object
        if (options.toolbar && options.toolbar.length) {
            var toolbarContainer = document.createElement('div');
            $(toolbarContainer).css('display', 'nonex');
            $(toolbarContainer).prop('class', 'jexcel-toolbar');

            $.each(options.toolbar, function(k, v) {
                if (v.type == 'i') {
                    if (v.method) {
                        var toolbarMethod = function() {
                            v.method($('#' + $.fn.jexcel.current), $.fn.jexcel('getSelectedCells'));
                        }
                    } else {
                        var toolbarMethod = function() {
                            $('#' + $.fn.jexcel.current).jexcel('setStyle', $.fn.jexcel('getSelectedCells'), v.k, v.v);
                        }
                    }
                    var toolbar = document.createElement('i');
                    $(toolbar).prop('class', 'jexcel-toolbar-item material-icons');
                    $(toolbar).on('click', toolbarMethod);
                    $(toolbar).html(v.content);
                    $(toolbarContainer).append(toolbar);
                } else if (v.type == 'select') {
                    var toolbar = document.createElement('select');
                    $(toolbar).prop('class', 'jexcel-toolbar-item');
                    $(toolbar).on('change', function() {
                        $('#' + $.fn.jexcel.current).jexcel('setStyle', $.fn.jexcel('getSelectedCells'), v.k, this.value);
                    });
                    if (typeof(v.v[0]) == 'string') {
                        $.each(v.v, function(k1, v1) {
                            var toolbarDropdownOption = document.createElement('option');
                            $(toolbarDropdownOption).prop('value', v1);
                            $(toolbarDropdownOption).html(v1);
                            $(toolbar).append(toolbarDropdownOption);
                        });
                    }
                    $(toolbarContainer).append(toolbar);
                } else if (v.type == 'spectrum') {
                    var toolbar = document.createElement('input');
                    $(toolbarContainer).append(toolbar);

                    $(toolbar).spectrum({
                        showButtons: false,
                        move: function(color) {
                            $('#' + $.fn.jexcel.current).jexcel('setStyle', $.fn.jexcel('getSelectedCells'), v.k, color.toHexString());
                        }
                    });
                }
            });
        }

        // Create header container
        var tableHeaderContainer = document.createElement('div');
        $(tableHeaderContainer).prop('class', 'jexcel-header');

        // Create content container
        var tableContentContainer = document.createElement('div');
        $(tableContentContainer).prop('class', 'jexcel-content');

        var tableHeader = document.createElement('table');
        $(tableHeader).prop('cellpadding', '0');
        $(tableHeader).prop('cellspacing', '0');
        $(tableHeader).prop('unselectable', 'yes');
        $(tableHeader).prop('onselectstart', 'return false');
        //$(table).prop('draggable', 'false');

        // Create header and body tags
        var thead = document.createElement('thead');
        var tbody = document.createElement('tbody');

        // Header
        $(thead).prop('class', 'jexcel_label');

        // Row
        var contentRow = '';
        var contentNested = '';
        var contentFilter = '';
        var contentWidth = '';

        // Create nested headers
        if (options.nestedHeaders && options.nestedHeaders.length > 0) {
            var createNestedHeader = function(nestedInformation) {
                // Create nested row
                contentRow = '<td></td>';

                var headerIndex = 0;
                for (var i = 0; i < nestedInformation.length; i++) {
                    // Default values
                    if (! nestedInformation[i].colspan) {
                        nestedInformation[i].colspan = 1;
                    }
                    if (! nestedInformation[i].align) {
                        nestedInformation[i].align = 'center';
                    }
                    if (! nestedInformation[i].title) {
                        nestedInformation[i].title = '';
                    }

                    // Classes container
                    var headerClasses = '';
                    // Header classes for this cell
                    for (var x = 0; x < nestedInformation[i].colspan; x++) {
                        if (headerClasses) {
                            headerClasses += ' ';
                        }
                        headerClasses += 'h' + headerIndex;
                        headerIndex++;
                    }
                    // Created the nested cell
                    contentRow += '<td class="'+ headerClasses +'" colspan="' + nestedInformation[i].colspan + '" align="' + nestedInformation[i].align + '">' + nestedInformation[i].title + '</td>';
                }

                return '<tr class="jexcel_nested">' + contentRow + '</tr>';
            }

            // Flexible way to handle nestedheaders
            if (options.nestedHeaders[0] && options.nestedHeaders[0][0]) {
                $.each(options.nestedHeaders, function(k, v) {
                    contentNested += createNestedHeader(v);
                });
            } else {
                contentNested += createNestedHeader(options.nestedHeaders);
            }
        }

        // Create headers
        contentRow = '<td></td>';

        // TODO: When the first or last column is hidden
        for (var i = 0; i < options.colHeaders.length; i++) {
            // Column properties
            var width = options.colWidths[i];
            var align = options.colAlignments[i] || 'center';
            var className = '';

            // Header class
            if (options.colHeaderClasses[i]) {
                className += options.colHeaderClasses[i];
            }

            // Display
            var display = (options.columns[i].type == 'hidden') ? ' style="display:none"' : '';

            // Default header cell properties
            var title = options.colHeaders[i];

            // Header
            var header = title || $.fn.jexcel('getColumnName', i);

            // Global column width
            contentWidth += '<col width="' + width + '">';

            // Create HTML row
            contentRow += '<td id="col-' + i + '" align="' + align +'" title="' + title + '" class="' + className + '"' + display + '>' + header + '</td>';

            // Filter columns
            contentFilter += '<td><input placeholder="Filter:"></td>';
        }

        // Populate header
        $(thead).append('<tr class="jexcel_headers_nested">' + contentNested + '</tr>');
        $(thead).append('<tr class="jexcel_headers">' + contentRow + '</tr>');
        //$(thead).append('<tr class="jexcel_filter"><td></td>' + contentFilter + '</tr>');

        // Append content
        $(tableHeader).append('<colgroup><col width="30">' + contentWidth + '</colgroup>');
        $(tableHeader).append(thead);
        $(tableHeaderContainer).append(tableHeader);

        // Create main table object
        var tableContent = document.createElement('table');
        $(tableContent).prop('cellpadding', '0');
        $(tableContent).prop('cellspacing', '0');
        $(tableContent).prop('unselectable', 'yes');
        $(tableContent).prop('onselectstart', 'return false');
        $(tableContent).append('<colgroup><col width="30">' + contentWidth + '</colgroup>');
        $(tableContent).append(tbody);
        $(tableContentContainer).html(tableContent);

        // Table overflow controls
        if ($.fn.jexcel.defaults[id].tableOverflow == true) {
            if ($.fn.jexcel.defaults[id].tableWidth) {
                $(tableContentContainer).css('overflow-x', 'scroll');
                $(tableContentContainer).css('max-width', $.fn.jexcel.defaults[id].tableWidth);
            }

            if ($.fn.jexcel.defaults[id].tableHeight) {
                $(tableContentContainer).css('overflow-y', 'scroll');
                $(tableContentContainer).css('max-height', $.fn.jexcel.defaults[id].tableHeight);

                // Adjust Header
                if ($.fn.jexcel.defaults[id].tableWidth) {
                    $(tableHeaderContainer).css('overflow-x', 'hidden');
                    $(tableHeaderContainer).css('max-width', parseInt($.fn.jexcel.defaults[id].tableWidth) - 17);
                }
            } else {
                if ($.fn.jexcel.defaults[id].tableWidth) {
                    $(tableHeaderContainer).css('overflow-x', 'hidden');
                    $(tableHeaderContainer).css('max-width', $.fn.jexcel.defaults[id].tableWidth);
                }
            }
        }

        // Main object
        $(this).html(toolbarContainer)
        $(this).append(tableHeaderContainer);
        $(this).append(tableContentContainer);

        // Add the corner square and textarea one time onlly
        if (! $('.jexcel_corner').length) {
            // Corner one for all sheets in a page
            var corner = document.createElement('div');
            $(corner).prop('class', 'jexcel_corner jexcel');
            $(corner).prop('id', 'jexcel_corner');

            // Disable copy
            if (options.selectionCopy == false) {
                $(corner).css('display', 'none');
            }

            // Hidden textarea copy and paste helper
            var textarea = document.createElement('textarea');
            $(textarea).prop('class', 'jexcel_textarea');
            $(textarea).prop('id', 'jexcel_textarea');

            // Contextmenu container
            var contextMenu = document.createElement('div');
            $(contextMenu).css('display', 'none');
            $(contextMenu).prop('class', 'jexcel_contextmenu');
            $(contextMenu).prop('id', 'jexcel_contextmenu');

            // Powered by
            var ads = document.createElement('div');
            $(ads).css('display', 'none');
            $(ads).prop('id', 'jexcel_about');
            $(ads).prop('class', 'jexcel_about');
            $(ads).html('<a href="https://github.com/paulhodel/jexcel">jExcel Spreadsheet</a>');

            // Append elements
            $('body').append(corner);
            $('body').append(textarea);
            $('body').append(contextMenu);
            $('body').append(ads);

            // Unselectable properties
            $(corner).prop('unselectable', 'on');
            $(corner).prop('onselectstart', 'return false');
            //$(corner).prop('draggable', 'false');

            // Prevent dragging on the corner object
            $.fn.jexcel.dragStartControls = function (e) {
                if ($(e.target).prop('id') == 'jexcel_corner') {
                   // return false;
                }
            }

            //$(document).on('dragstart', $.fn.jexcel.dragStartControls);

            // Corner persistence and other helpers
            $.fn.jexcel.selectedCorner = false;
            $.fn.jexcel.selectedHeader = null;
            $.fn.jexcel.resizeColumn = null;

            // Context menu
            $.fn.jexcel.contextMenuControls = function (e) {
                if ($($.fn.jexcel.selectedCell).hasClass('edition')) {
                    e.preventDefault();
                } else {
                    // Hide jExcel context menu if is open
                    if ($("#jexcel_contextmenu").css('display') == 'block') {
                        $("#jexcel_contextmenu").css('display', 'none')
                    }

                    if ($.fn.jexcel.current) {
                        // Check if the click was in an jexcel element
                        var isJexcel = $(e.target).parents('.jexcel').length ? true : false;

                        // Table found
                        if (isJexcel) {
                            // The id is depending on header and body
                            if ($(e.target).parent().is('.jexcel_headers')) {
                                var o = $(e.target).prop('id');
                                var h = true;
                            } else {
                                var o = $(e.target).parent().prop('id');
                                var h = false;
                            }

                            if (o) {
                                o = o.split('-');
                                contextMenuContent = '';
                                // Custom context menu
                                if (typeof($.fn.jexcel.defaults[$.fn.jexcel.current].contextMenu) == 'function') {
                                    contextMenuContent = $.fn.jexcel.defaults[$.fn.jexcel.current].contextMenu(o[0], o[1], e);
                                } else {
                                    // Default context menu for the columns
                                    contextMenuContent = $.fn.jexcel('contextMenu', o[0], o[1], e);
                                }

                                // Show context menu
                                if (contextMenuContent) {
                                    // Contextmenu content
                                    $("#jexcel_contextmenu").html(contextMenuContent);

                                    // Show jexcel context menu
                                    $("#jexcel_contextmenu").css({
                                        display: 'block',
                                        top: e.pageY + 'px',
                                        left: e.pageX + 'px'
                                    });

                                    // Avoid the real one
                                    e.preventDefault();
                                }
                            }
                        }
                    }
                }
            }

            $(document).on("contextmenu", $.fn.jexcel.contextMenuControls);

            // Mouse wheel
            $.fn.jexcel.mouseWheelControls = function (e) {
                // Hide context menu
                $(".jexcel_contextmenu").css('display', 'none');
            }

            $(document).on('mousewheel', $.fn.jexcel.mouseWheelControls);

            // Global mouse click down controles
            $.fn.jexcel.mouseDownControls = function (e) {
                if (e.which != 3) {
                    // Click on corner icon
                    if (e.target.id == 'jexcel_corner') {
                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                            $.fn.jexcel.selectedCorner = true;
                        }
                    } else {
                        // Check if the click was in an jexcel element
                        var jexcelTable = $(e.target).parents('.jexcel');

                        // Table found
                        if ($(jexcelTable).length > 0) {
                            // Get id
                            var current = $(jexcelTable).prop('id');

                            // Remove selection from any other jexcel if applicable
                            if ($.fn.jexcel.current) {
                                if ($.fn.jexcel.current != current) {
                                    $('#' + $.fn.jexcel.current).jexcel('updateSelection');
                                }
                            }

                            // Mark as current
                            $.fn.jexcel.current = current;

                            // Header found
                            if ($(e.target).parent().is('.jexcel_headers')) {
                                var o = $(e.target).prop('id');
                                if (o) {
                                    o = o.split('-');

                                    if ($.fn.jexcel.selectedHeader && (e.shiftKey || e.ctrlKey)) {
                                        var d = $($.fn.jexcel.selectedHeader).prop('id').split('-');
                                    } else {
                                        // Update selection single column
                                        var d = $(e.target).prop('id').split('-');
                                        // Keep track of which header was selected first
                                        $.fn.jexcel.selectedHeader = $(e.target);
                                        $.fn.jexcel.selectedRow = null;
                                    }

                                    // Update cursor
                                    if ($(e.target).outerWidth() - e.offsetX < 8) {
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].columnResize == true) {
                                            // Resize helper
                                            $.fn.jexcel.resizeColumn = {
                                                mousePosition: e.pageX,
                                                column:o[1],
                                                width:parseInt($(e.target).css('width')),
                                            }
                                            // Border indication
                                            $('#' + $.fn.jexcel.current + ' .c' + o[1]).addClass('resizing');
                                            $('#' + $.fn.jexcel.current + ' .h' + o[1]).addClass('resizing');
                                            $('#' + $.fn.jexcel.current + ' #col-' + o[1]).addClass('resizing');

                                            // Remove selected cells
                                            $('#' + $.fn.jexcel.current).jexcel('updateSelection');
                                            $('.jexcel_corner').css('left', '-200px');
                                        }
                                    } else {
                                        // Get cell objects
                                        var o1 = $('#' + $.fn.jexcel.current).find('#' + o[1] + '-0');
                                        var o2 = $('#' + $.fn.jexcel.current).find('#' + d[1] + '-' + parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].data.length - 1));

                                        // Update selection
                                        $('#' + $.fn.jexcel.current).jexcel('updateSelection', o1, o2, 1);

                                        // Selected cell will be the first in the row
                                        $.fn.jexcel.selectedCell = $(o1);
                                    }
                                }
                            } else {
                                $.fn.jexcel.selectedHeader = false;
                            }

                            // Nested headers
                            if ($(e.target).parent().is('.jexcel_nested')) {
                                // Reference class headers
                                var ob = $(e.target).prop('class').split(' ');
                                var st = '';

                                // Get the columns should be included in the search
                                $.each(ob, function(k, v) {
                                    if (v.substr(0,1) == 'h') {
                                        if (st) {
                                            st += ', ';
                                        }
                                        st += '.c' + v.substr(1);
                                    }
                                });

                                if (st) {
                                    // Search all cells in the defined search
                                    var ob = $('#' + $.fn.jexcel.current).find(st);
                                    var o1 = $(ob).first();
                                    var o2 = $(ob).last();

                                    // Update selection
                                    $('#' + $.fn.jexcel.current).jexcel('updateSelection', o1, o2, 1);

                                    // Selected cell will be the first in the row
                                    $.fn.jexcel.selectedCell = $(o1);
                                }
                            }

                            // Body found
                            if ($(e.target).parents('.jexcel-content').length) {
                                var target = (e.target.tagName == 'TD') ? e.target : $(e.target).parents('td');
                                if ($(target).length && $(target).prop('id')) {
                                    // Update row label selection
                                    if ($(e.target).is('.jexcel_label')) {
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].rowDrag == true && $(e.target).outerWidth() - e.offsetX < 8) {
                                            // Reset selection
                                            $('#' + $.fn.jexcel.current).jexcel('resetSelection');
                                            // Mark which row we are dragging
                                            $.fn.jexcel.dragRowFrom = $(e.target).parent().prop('id');
                                            $.fn.jexcel.dragRowOver = $(e.target).parent().prop('id');
                                            // Visual row we are dragging
                                            $(e.target).parent().find('td').css('background-color', 'rgba(0,0,0,0.1)');
                                        } else {
                                            var o = $(e.target).parent().prop('id').split('-');

                                            if ($.fn.jexcel.selectedRow && (e.shiftKey || e.ctrlKey)) {
                                                // Updade selection multi columns
                                                var d = $($.fn.jexcel.selectedRow).prop('id').split('-');
                                            } else {
                                                // Update selection single column
                                                var d = $(e.target).parent().prop('id').split('-');
                                                // Keep track of which header was selected first
                                                $.fn.jexcel.selectedRow = $(e.target).parent();
                                                $.fn.jexcel.selectedHeader = null;
                                            }

                                            // Get cell objects
                                            var o1 = $('#' + $.fn.jexcel.current).find('#0-' + o[1]);
                                            var o2 = $('#' + $.fn.jexcel.current).find('#' + parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].columns.length - 1) + '-' + d[1]);

                                            // Update selection
                                            $('#' + $.fn.jexcel.current).jexcel('updateSelection', o1, o2);

                                            // Selected cell will be the first in the row
                                            $.fn.jexcel.selectedCell = $(o2);
                                        }
                                    } else {
                                        // Update cell selection
                                        if (! $(e.target).parents('.edition').length) {
                                            if (! $.fn.jexcel.selectedCell || ! e.shiftKey) {
                                                $.fn.jexcel.selectedCell = $(target);
                                            }
                                            $('#' + $.fn.jexcel.current).jexcel('updateSelection', $.fn.jexcel.selectedCell, $(target));
                                        } else {
                                            if ($(target) != $.fn.jexcel.selectedCell) {
                                                $.fn.jexcel.selectedCell = $(target);
                                                $('#' + $.fn.jexcel.current).jexcel('updateSelection', $.fn.jexcel.selectedCell, $(target));
                                            }
                                        }

                                        // Edition 
                                        $.fn.jexcel.touchControls = setTimeout(function() {
                                            if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                                if ($(e.target).parents('.jexcel-content').length) {
                                                    var target = (e.target.tagName == 'TD') ? e.target : $(e.target).parents('td');

                                                    // Open editor action
                                                    if ($(target).is('.highlight')) {
                                                        $('#' + $.fn.jexcel.current).jexcel('openEditor', $(target), null, e);
                                                    }
                                                }
                                            }
                                        }, 500);

                                        // No full row selected
                                        $.fn.jexcel.selectedRow = null;
                                    }
                                }
                            }
                        } else {
                            // Check if the object is in the jexcel domain
                            if (! $(e.target).parents('.jexcel, .jexcel_contextmenu').length) {
                                // Keep selection if main scrollbar is selected
                                if (e.target != $('html').get(0)) {
                                    $('#' + $.fn.jexcel.current).jexcel('resetSelection');
                                }
                            }
                        }
                    }
                }
            }

            $(document).on('mousedown touchstart', $.fn.jexcel.mouseDownControls);

            // Global mouse click up controles
            $.fn.jexcel.mouseUpControls = function (e) {
                // The click in array means open the dropdown
                if (e.target.id == 'jexcel_arrow') {
                    // If not current table is selected
                    if (! $.fn.jexcel.current) {
                        $.fn.jexcel.current = $(e.target).parents('.jexcel').parent().prop('id');
                    }
                    // Mark this cell is the current cell cursor
                    $.fn.jexcel.selectedCell = $(e.target).parent().parent();
                    // Update selection for this cell
                    $('#' + $.fn.jexcel.current).jexcel('updateSelection', $.fn.jexcel.selectedCell, $.fn.jexcel.selectedCell);
                    // Open editor
                    $('#' + $.fn.jexcel.current).jexcel('openEditor', $.fn.jexcel.selectedCell, null, e);
                } else {
                    if (e.which != 3) {
                        // Hide context menu
                        $("#jexcel_contextmenu").css('display', 'none');
                    }

                    // Cancel any corner selection
                    $.fn.jexcel.selectedCorner = false;

                    // Update cell size
                    if ($.fn.jexcel.resizeColumn) {
                        // Columns to be updated
                        $('#' + $.fn.jexcel.current + ' .h' + $.fn.jexcel.resizeColumn.column).removeClass('resizing');
                        $('#' + $.fn.jexcel.current + ' .c' + $.fn.jexcel.resizeColumn.column).removeClass('resizing');
                        $('#' + $.fn.jexcel.current + ' #col-' + $.fn.jexcel.resizeColumn.column).removeClass('resizing');

                        // New width
                        var newWidth = $('#' + $.fn.jexcel.current + ' .jexcel-header col').eq(parseInt($.fn.jexcel.resizeColumn.column) + 1).prop('width');

                        // Update width
                        $('#' + $.fn.jexcel.current).jexcel('setWidth', $.fn.jexcel.resizeColumn.column, newWidth);

                        // On resize
                        if (typeof($.fn.jexcel.defaults[$.fn.jexcel.current].onresize) == 'function') {
                            $.fn.jexcel.defaults[$.fn.jexcel.current].onresize($(this), $.fn.jexcel.resizeColumn.column, newWidth, $.fn.jexcel.resizeColumn.width);
                        }

                        // Reset resizing helper
                        $.fn.jexcel.resizeColumn = null;
                    }

                    // Data to be copied
                    var selection = $('#' + $.fn.jexcel.current).find('tbody td.selection');

                    if ($(selection).length > 0) {
                        // First and last cells
                        var o = $(selection[0]).prop('id').split('-');
                        var d = $(selection[selection.length - 1]).prop('id').split('-');

                        // Copy data
                        $('#' + $.fn.jexcel.current).jexcel('copyData', o, d);

                        // Remove selection
                        $(selection).removeClass('selection selection-left selection-right selection-top selection-bottom');
                    }
                }

                // Execute the final move TODO - finish this...
                if ($.fn.jexcel.dragRowFrom) {
                    if ($.fn.jexcel.dragRowFrom != $.fn.jexcel.dragRowOver) {
                        var fr = $.fn.jexcel.dragRowFrom.split('-');
                        var to = $.fn.jexcel.dragRowOver.split('-');
                        $('#' + $.fn.jexcel.current).jexcel('moveRow', fr[1], to[1]);
                    }

                    // Remove style
                    $('#' + $.fn.jexcel.dragRowFrom).css('cursor', '');
                    $('#' + $.fn.jexcel.dragRowOver).css('cursor', '');
                    $('#' + $.fn.jexcel.current + ' #' + $.fn.jexcel.dragRowOver).parent().find('td').css('background-color', '');
                }

                $.fn.jexcel.dragRowFrom = null;
                $.fn.jexcel.dragRowOver = null;

                // Cancel any timing
                if ($.fn.jexcel.touchControls) {
                    clearTimeout($.fn.jexcel.touchControls);
                }
            }

            $(document).on('mouseup', $.fn.jexcel.mouseUpControls);

            // Touch control
            $.fn.jexcel.touchControls = null

            $.fn.jexcel.touchCancelControls  = function() {
                if ($.fn.jexcel.touchControls) {
                    clearTimeout($.fn.jexcel.touchControls);
                }
            }

            $(document).on('touchend touchcancel touchmove', $.fn.jexcel.touchCancelControls);

            // Double click
            $.fn.jexcel.doubleClickControls = function (e) {
                // Jexcel is selected
                if ($.fn.jexcel.current) {
                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                        // Corner action
                        if (e.target.id == 'jexcel_corner') {
                            var selection = $('#' + $.fn.jexcel.current).find('tbody td.highlight');
                            // Any selected cells
                            if (typeof(selection) == 'object') {
                                // Get selected cells
                                var o = $(selection[0]).prop('id').split('-');
                                var d = $(selection[selection.length - 1]).prop('id').split('-');
                                // Double click copy
                                o[1] = parseInt(d[1]) + 1;
                                d[1] = parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].data.length);
                                // Do copy
                                $('#' + $.fn.jexcel.current).jexcel('copyData', o, d);
                            }
                        } else if ($(e.target).parents('.jexcel-content').length) {
                            var target = (e.target.tagName == 'TD') ? e.target : $(e.target).parents('td');

                            // Open editor action
                            if ($(target).is('.highlight')) {
                                $('#' + $.fn.jexcel.current).jexcel('openEditor', $(target), null, e);
                            }
                        }
                    }

                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].columnSorting == true) {
                        // Header found
                        if ($(e.target).parent().parent().is('thead')) {
                            var o = $(e.target).prop('id');
                            if (o) {
                                o = $(e.target).prop('id').split('-');
                                // Update order
                                $('#' + $.fn.jexcel.current).jexcel('orderBy', o[1]);
                            }
                        }
                    }
                }
            }

            $(document).on('dblclick', $.fn.jexcel.doubleClickControls);

            // Mouse move controls
            $.fn.jexcel.mouseMoveControls = function (e) {
                if ($.fn.jexcel.current) {
                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].columnResize == true) {
                        // Resizing is ongoing
                        if ($.fn.jexcel.resizeColumn) {
                           var width = e.pageX - $.fn.jexcel.resizeColumn.mousePosition;

                           if ($.fn.jexcel.resizeColumn.width + width > 0) {
                               var tempWidth = $.fn.jexcel.resizeColumn.width + width;
                               $('#' + $.fn.jexcel.current).find('.jexcel-header col').eq(parseInt($.fn.jexcel.resizeColumn.column) + 1).prop('width', tempWidth);
                           }
                        } else {
                            // Header found
                            if ($(e.target).parent().is('.jexcel_headers')) {
                                // Update cursor
                                if ($(e.target).outerWidth() - e.offsetX < 8 && $(e.target).prop('id') != '') {
                                    $(e.target).css('cursor', 'col-resize');
                                } else {
                                    if ($(e.target).css('cursor') == 'col-resize') {
                                        $(e.target).css('cursor', '');
                                    }
                                }
                            }
                        }
                    }

                    // Body found
                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].rowDrag == true) {
                        if ($(e.target).parent().parent().is('tbody')) {
                            // Update row label selection
                            if ($(e.target).is('.jexcel_label')) {
                                if ($(e.target).outerWidth() - e.offsetX < 8) {
                                    $(e.target).css('cursor', 'all-scroll');
                                } else {
                                    if ($(e.target).css('cursor') == 'all-scroll') {
                                        $(e.target).css('cursor', '');
                                    }
                                }
                            }
                        }
                    }
                }

                // Keeping visual indication
                if ($.fn.jexcel.dragRowFrom) {
                    $('body').css('cursor', 'all-scroll');
                } else {
                    if ($('body').css('cursor') == 'all-scroll') {
                        $('body').css('cursor', '');
                    }
                }

                // Cancel timing for edition
                if ($.fn.jexcel.touchControls) {
                    clearTimeout($.fn.jexcel.touchControls);
                }
            }

            $(document).on('mousemove', $.fn.jexcel.mouseMoveControls);

            // Mouse over controls
            $.fn.jexcel.mouseOverControls = function (e) {
                // No resizing is ongoing
                if (! $.fn.jexcel.resizeColumn) {
                    // Get jexcel table
                    var table = $(e.target).closest('.jexcel');

                    // If the user is in the current table
                    if ($.fn.jexcel.current == $(table).prop('id')) {
                        // Header found
                        if ($(e.target).parent().is('.jexcel_headers')) {
                            if ($.fn.jexcel.selectedHeader) {
                                // Updade selection
                                if (e.buttons) {
                                    var o = $($.fn.jexcel.selectedHeader).prop('id');
                                    var d = $(e.target).prop('id');
                                    if (o && d) {
                                        o = o.split('-');
                                        d = d.split('-');
                                        // Get cell objects
                                        var o1 = $('#' + $.fn.jexcel.current).find('#' + o[1] + '-0');
                                        var o2 = $('#' + $.fn.jexcel.current).find('#' + d[1] + '-' + parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].data.length - 1));
                                        // Update selection
                                        $('#' + $.fn.jexcel.current).jexcel('updateSelection', o1, o2);
                                    }
                                }
                            }
                        }

                        // Body found
                        if ($(e.target).parent().parent().is('tbody')) {
                            // Update row label selection
                            if ($(e.target).is('.jexcel_label')) {
                                if ($.fn.jexcel.selectedRow) {
                                    // Updade selection
                                    if (e.buttons) {
                                        var o = $($.fn.jexcel.selectedRow).prop('id');
                                        var d = $(e.target).parent().prop('id');
                                        if (o && d) {
                                            o = o.split('-');
                                            d = d.split('-');
                                            // Get cell objects
                                            var o1 = $('#' + $.fn.jexcel.current).find('#0-' + o[1]);
                                            var o2 = $('#' + $.fn.jexcel.current).find('#' + parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].columns.length - 1) + '-'  + d[1]);
                                            // Update selection
                                            $('#' + $.fn.jexcel.current).jexcel('updateSelection', o1, o2);
                                        }
                                    }
                                } else if ($.fn.jexcel.dragRowFrom) {
                                    // Remove previous row visual background
                                    $('#' + $.fn.jexcel.current + ' #' + $.fn.jexcel.dragRowOver).parent().find('td').css('background-color', '');
                                    // Add new row visual background
                                    $(e.target).parent().find('td').css('background-color', 'rgba(0,0,0,0.1)');
                                    // Keep over reference
                                    $.fn.jexcel.dragRowOver = $(e.target).parent().prop('id');
                                }
                            } else {
                                if ($.fn.jexcel.selectedCell) {
                                    if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                        if ($.fn.jexcel.selectedCorner == true) {
                                            // Copy option
                                            $('#' + $.fn.jexcel.current).jexcel('updateCornerSelection', $(e.target));
                                        } else {
                                            // Updade selection
                                            if (e.buttons) {
                                                $('#' + $.fn.jexcel.current).jexcel('updateSelection', $.fn.jexcel.selectedCell, $(e.target));
                                            }
                                        }
                                    }
                                }

                                // On mouse over event for cells
                                if ($.fn.jexcel.ignoreEvents != true) {
                                    if ($.fn.jexcel.defaults[id].onmouseover) {
                                        if (typeof($.fn.jexcel.defaults[id].onmouseover) == 'function') {
                                            $.fn.jexcel.defaults[id].onmouseover($('#' + $.fn.jexcel.current), $(e.target));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            $(document).on('mouseover', $.fn.jexcel.mouseOverControls);

            // Paste Controls - IE Compatibility
            $.fn.jexcel.pasteControls = function (e) {
                if ($.fn.jexcel.selectedCell) {
                    if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                        if (e.originalEvent) {
                            if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                $('#' + $.fn.jexcel.current).jexcel('paste', $.fn.jexcel.selectedCell, e.originalEvent.clipboardData.getData('text'));
                            }
                            e.preventDefault();
                        }
                    }
                }
            }

            $(document).on('paste', $.fn.jexcel.pasteControls);

            // Key down controls
            $.fn.jexcel.keyDownControls = function(e) {
                if ($.fn.jexcel.current) {
                    // Support variables
                    var cell = null;
                    // Get current cell
                    if ($.fn.jexcel.selectedCell) {
                        columnId = $($.fn.jexcel.selectedCell).prop('id').split('-');

                        // Which key
                        if (e.which == 37) {
                            // Left arrow
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if (e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().find('td').not('.jexcel_label').first();
                                } else {
                                    if ($.fn.jexcel('col') > 0) {
                                        cell = $($.fn.jexcel.selectedCell).prev();
                                    }
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 39) {
                            // Right arrow
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if (e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().find('td').last();
                                } else {
                                    cell = $($.fn.jexcel.selectedCell).next();
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 38) {
                            // Top arrow
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if (e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().parent().find('tr').first().find('#' + columnId[0] + '-' + 0);
                                } else {
                                    cell = $($.fn.jexcel.selectedCell).parent().prev().find('#' + columnId[0] + '-' + (columnId[1] - 1));
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 40) {
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                // Bottom arrow
                                if (e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().parent().find('tr').last().find('#' + columnId[0] + '-' + ($.fn.jexcel.defaults[$.fn.jexcel.current].data.length - 1));
                                } else {
                                    cell = $($.fn.jexcel.selectedCell).parent().next().find('#' + columnId[0] + '-' + (parseInt(columnId[1]) + 1));
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 36) {
                            // Home
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if (e.shiftKey && e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().parent().find('tr').first().find('td').not('.jexcel_label').first();
                                } else {
                                    cell = $($.fn.jexcel.selectedCell).parent().find('td').not('.jexcel_label').first();
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 35) {
                            // End
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if (e.shiftKey && e.ctrlKey) {
                                    cell = $($.fn.jexcel.selectedCell).parent().parent().find('tr').last().find('td').last();
                                } else {
                                    cell = $($.fn.jexcel.selectedCell).parent().find('td').last();
                                }
                                e.preventDefault();
                            }
                        } else if (e.which == 27) {
                            // Escape
                            if ($($.fn.jexcel.selectedCell).hasClass('edition')) {
                                // Exit without saving
                                $('#' + $.fn.jexcel.current).jexcel('closeEditor', $($.fn.jexcel.selectedCell), false);
                            }
                        } else if (e.which == 13) {
                            // Enter
                            if ($($.fn.jexcel.selectedCell).hasClass('edition')) {
                                // Exit saving data
                                if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'calendar') {
                                    $('#' + $.fn.jexcel.current).find('editor').jcalendar('close', 1);
                                } else if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'dropdown' || $.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'autocomplete') {
                                    // Do nothing
                                } else {
                                    // Alt enter -> do not close editor
                                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].wordWrap == true && e.altKey) {
                                        // Add new line to the editor
                                        var editorTextarea = $($.fn.jexcel.selectedCell).find('.editor');
                                        var editorValue = $(editorTextarea).val();
                                        var editorIndexOf = $(editorTextarea).prop('selectionStart');
                                        editorValue = editorValue.slice(0, editorIndexOf) + "\n" + editorValue.slice(editorIndexOf);
                                        $(editorTextarea).val(editorValue);
                                        $(editorTextarea).focus();
                                        $(editorTextarea).prop('selectionStart', editorIndexOf + 1);
                                        $(editorTextarea).prop('selectionEnd', editorIndexOf + 1);
                                    } else {
                                        $('#' + $.fn.jexcel.current).jexcel('closeEditor', $($.fn.jexcel.selectedCell), true);
                                    }
                                }
                            }

                            // If not edition check if the selected cell is in the last row
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowInsertRow == true) {
                                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowManualInsertRow == true) {
                                        if (columnId[1] == $.fn.jexcel.defaults[$.fn.jexcel.current].data.length - 1) {
                                            // New record in case selectedCell in the last row
                                            $('#' + $.fn.jexcel.current).jexcel('insertRow');
                                        }
                                    }
                                }
                                // Go to the next line
                                cell = $($.fn.jexcel.selectedCell).parent().next().find('#' + columnId[0] + '-' + (parseInt(columnId[1]) + 1));
                                e.preventDefault();
                            }
                        } else if (e.which == 9) {
                            // Tab
                            if ($($.fn.jexcel.selectedCell).hasClass('edition')) {
                                // Exit saving data
                                if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'calendar') {
                                    $('#' + $.fn.jexcel.current).find('editor').jcalendar('close', 1)
                                } else {
                                    $('#' + $.fn.jexcel.current).jexcel('closeEditor', $($.fn.jexcel.selectedCell), true);
                                }
                            }
                            // Tab key - Get the id of the selected cell
                            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowInsertColumn == true) {
                                if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowManualInsertColumn == true) {
                                    if (columnId[0] == $.fn.jexcel.defaults[$.fn.jexcel.current].data[0].length - 1) {
                                        // New record in case selectedCell in the last column
                                        $('#' + $.fn.jexcel.current).jexcel('insertColumn');
                                    }
                                }
                            }
                            // Highlight new column
                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                               if (e.shiftKey) {
                                   if ($.fn.jexcel('col') > 0) {
                                       cell = $($.fn.jexcel.selectedCell).prev();
                                   }
                               } else {
                                   cell = $($.fn.jexcel.selectedCell).next();
                               }
                            }
                            e.preventDefault();
                        } else if (e.which == 46) {
                            // Delete (erase cell in case no edition is running)
                            if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                    if ($.fn.jexcel.selectedRow) {
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowDeleteRow == true) {
                                            if (confirm('Are you sure to delete the selected rows?')) {
                                                 $('#' + $.fn.jexcel.current).jexcel('deleteRow');
                                                 $('#' + $.fn.jexcel.current).jexcel('resetSelection');
                                            }
                                        }
                                    } else if ($.fn.jexcel.selectedHeader) {
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowDeleteColumn == true) {
                                            if (confirm('Are you sure to delete the selected columns?')) {
                                                var columns = $('#' + $.fn.jexcel.current).find('thead.jexcel_label').find('.selected');
                                                $('#' + $.fn.jexcel.current).jexcel('deleteColumn', parseInt($(columns[0]).prop('id').split('-')[1]), parseInt(columns.length));
                                            }
                                        }
                                    } else {
                                        // Change value
                                        $('#' + $.fn.jexcel.current).jexcel('setValue', $('#' + $.fn.jexcel.current).find('.highlight'), '');
                                    }
                                }
                            }
                        } else {
                            if (e.metaKey && ! e.shiftKey && ! e.ctrlKey) {
                                if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                    if (e.which == 67) {
                                        // Command + C, Mac
                                        $('#' + $.fn.jexcel.current).jexcel('copy', true);
                                        e.preventDefault();
                                    }
                                }
                            } else if (! e.ctrlKey) {
                                if ($.fn.jexcel.selectedCell) {
                                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                        // If is not readonly
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type != 'readonly') {
                                            // Start edition in case a valid character.
                                            if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                                // Characters able to start a edition
                                                if (e.keyCode == 32) {
                                                    // Space
                                                    if ($.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'checkbox' || $.fn.jexcel.defaults[$.fn.jexcel.current].columns[columnId[0]].type == 'radio') {
                                                        var checkboxCurrentVal = $('#' + $.fn.jexcel.current).jexcel('getValue', $.fn.jexcel.selectedCell);
                                                        $('#' + $.fn.jexcel.current).jexcel('setValue', $.fn.jexcel.selectedCell, checkboxCurrentVal == 1 ? false : true);
                                                        e.preventDefault();
                                                    } else {
                                                        $('#' + $.fn.jexcel.current).jexcel('openEditor', $($.fn.jexcel.selectedCell), true, e);
                                                    }
                                                } else if (e.keyCode == 110 || (e.keyCode >= 48 && e.keyCode <= 57) || (e.keyCode >= 65 && e.keyCode <= 90) || (e.keyCode >= 96 && e.keyCode <= 105) || (e.keyCode >= 186 && e.keyCode <= 190)) {
                                                    $('#' + $.fn.jexcel.current).jexcel('openEditor', $($.fn.jexcel.selectedCell), true, e);
                                                } else if (e.keyCode == 113) {
                                                    $('#' + $.fn.jexcel.current).jexcel('openEditor', $($.fn.jexcel.selectedCell), false, e);
                                                    e.preventDefault();
                                                }
                                            }
                                        }
                                    }
                                }
                            } else if (! e.shiftKey && e.ctrlKey) {
                                if (! $($.fn.jexcel.selectedCell).hasClass('edition')) {
                                    if (e.which == 65) {
                                        // Ctrl + A
                                        var o1 = parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].columns.length - 1);
                                        var o2 = parseInt($.fn.jexcel.defaults[$.fn.jexcel.current].data.length - 1);
                                        var o = $('#' + $.fn.jexcel.current).find('#0-0');
                                        var t = $('#' + $.fn.jexcel.current).find('#' + o1 + '-' + o2);
                                        $('#' + $.fn.jexcel.current).jexcel('updateSelection', o, t);
                                        // Prevent page selection
                                        e.preventDefault();
                                    } else if (e.which == 83) {
                                        // Ctrl + S
                                        $('#' + $.fn.jexcel.current).jexcel('download');
                                        // Prevent page selection
                                        e.preventDefault();
                                    } else if (e.which == 89) {
                                        // Ctrl + Y
                                        if (!$($.fn.jexcel.selectedCell).hasClass('edition')) {
                                            $('#' + $.fn.jexcel.current).jexcel('redo');
                                        }
                                        e.preventDefault();
                                    } else if (e.which == 90) {
                                        // Ctrl + Z
                                        if (!$($.fn.jexcel.selectedCell).hasClass('edition')) {
                                            $('#' + $.fn.jexcel.current).jexcel('undo');
                                        }
                                        e.preventDefault();
                                    } else if (e.which == 67) {
                                        // Ctrl + C
                                        $('#' + $.fn.jexcel.current).jexcel('copy', true);
                                        e.preventDefault();
                                    } else if (e.which == 88) {
                                        // Ctrl + X
                                        if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                            $('#' + $.fn.jexcel.current).jexcel('cut');
                                        } else {
                                            $('#' + $.fn.jexcel.current).jexcel('copy', true);
                                        }
                                        e.preventDefault();
                                    } else if (e.which == 86) {
                                        // Ctrl + V
                                        if (window.clipboardData) {
                                            if ($.fn.jexcel.defaults[$.fn.jexcel.current].editable == true) {
                                                $('#' + $.fn.jexcel.current).jexcel('paste', $.fn.jexcel.selectedCell, window.clipboardData.getData('Text'));
                                            }
                                            e.preventDefault();
                                        }
                                    }
                                }
                            }
                        }

                        // Arrows control
                        if (cell) {
                            // Control selected cell
                            if ($(cell).length > 0) {
                                // In case of a multiple cell selection
                                if (e.shiftKey && e.which != 9) {
                                    // Keep first selected cell
                                    if (! $.fn.jexcel.referenceCell) {
                                        $.fn.jexcel.referenceCell = $.fn.jexcel.selectedCell;
                                    }

                                    // Origin cell
                                    var o = $.fn.jexcel.referenceCell;
                                } else {
                                    // Remove previous cell
                                    $.fn.jexcel.referenceCell = null;

                                    // Origin cell
                                    var o = cell;
                                }

                                // Target cell
                                var d = cell;

                                // Current cell
                                $.fn.jexcel.selectedCell = d;

                                // Update selection
                                $('#' + $.fn.jexcel.current).jexcel('updateSelection', o, d);
                            }
                        }
                    }
                }
            }

            $(document).on('keydown', $.fn.jexcel.keyDownControls);
        }

        $(this).children('.jexcel-content').on('scroll', function() {
            $(this).jexcel('updateCornerPosition');
            $(tableHeaderContainer).scrollLeft($(this).scrollLeft());
        });

        // Load data
        $(this).jexcel('setData', $.fn.jexcel.defaults[id].data);

        // Update values
        var ignoreEvents = $.fn.jexcel.ignoreEvents ? true : false;
        var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

        // Disabled events and history
        $.fn.jexcel.ignoreEvents = true;
        $.fn.jexcel.ignoreHistory = true;

        // Load style
        if ($.fn.jexcel.defaults[id].style) {
            $(this).jexcel('setStyle', $.fn.jexcel.defaults[id].style);
        }

        // Load meta
        if ($.fn.jexcel.defaults[id].meta) {
            $(this).jexcel('setMeta', $.fn.jexcel.defaults[id].meta);
        }

        // Restore events and history flag
        $.fn.jexcel.ignoreEvents = ignoreEvents;
        $.fn.jexcel.ignoreHistory = ignoreHistory;
    },

    /**
     * Set data
     * 
     * @param array data In case no data is sent, default is reloaded
     * @return void
     */
    setData : function(data, ignoreSpare) {
        // Id
        var id = $(this).prop('id');

        // Update data
        if (data) {
            if (typeof(data) == 'string') {
                data = JSON.parse(data);
            }

            $.fn.jexcel.defaults[id].data = data;
        }

        // Create history track array
        $.fn.jexcel.defaults[id].history = [];
        $.fn.jexcel.defaults[id].historyIndex = -1;

        // Adjust minimal dimensions
        var size_i = $.fn.jexcel.defaults[id].colHeaders.length;
        var size_j = $.fn.jexcel.defaults[id].data.length;
        var min_i = $.fn.jexcel.defaults[id].minDimensions[0];
        var min_j = $.fn.jexcel.defaults[id].minDimensions[1];
        var max_i = min_i > size_i ? min_i : size_i;
        var max_j = min_j > size_j ? min_j : size_j;

        for (j = 0; j < max_j; j++) {
            for (i = 0; i < max_i; i++) {
                if ($.fn.jexcel.defaults[id].data[j] == undefined) {
                    $.fn.jexcel.defaults[id].data[j] = [];
                }

                if ($.fn.jexcel.defaults[id].data[j][i] == undefined) {
                    $.fn.jexcel.defaults[id].data[j][i] = '';
                }
            }
        }

        // Dynamic columns
        $.fn.jexcel.defaults[id].formula = [];

        // Values
        $.fn.jexcel.defaults[id].values = [];

        // Data container
        var tbody = $(this).find('tbody');

        // Reset data
        $(tbody).html('');

        // Records
        var records = [];

        // Create cells
        for (j = 0; j < $.fn.jexcel.defaults[id].data.length; j++) {
            // New line of data to be append in the table
            tr = document.createElement('tr');
            $(tr).prop('id', 'row-' + parseInt(j));

            // Index column
            $(tr).append('<td class="jexcel_label">' + parseInt(j + 1) + '</td>'); 

            // Data columns
            for (i = 0; i < $.fn.jexcel.defaults[id].colHeaders.length; i++) {
                // New column of data to be append in the line
                td = $(this).jexcel('createCell', i, j);

                // Add column to the row
                $(tr).append(td);

                // Cell data
                records.push({
                    col:i,
                    row:j,
                    cell: $(td),
                    newValue: $.fn.jexcel.defaults[id].data[j][i],
                    oldValue: '',
                });
            }

            // Add row to the table body
            $(tbody).append(tr);
        }

        // Update all records
        $(this).jexcel('loadCells', records, true);

        // Update all cells with formulas
        $(this).jexcel('updateAllCellsWithFormulas');

        // New data available
        if (data) {
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof($.fn.jexcel.defaults[id].onload) == 'function') {
                    $.fn.jexcel.defaults[id].onload($(this));
                }
            }
        }
    },

    /**
     * Update table settings helper. Update cells after loading
     * 
     * @param methods
     * @return void
     */
    updateSettings : function(options) {
        var main = $(this);

        // Id
        var id = $(this).prop('id');

        // Keep options
        if (options) {
            $.fn.jexcel.defaults[id].updateSettingsOptions = options;
        } else {
            if ($.fn.jexcel.defaults[id].updateSettingsOptions) {
                options = $.fn.jexcel.defaults[id].updateSettingsOptions;
            }
        }

        // Go through all cells
        if (options) {
            // Get all cells form 
            var cells = $(this).find('tbody td').not('.jexcel_label');
            // Existing methods
            if (typeof(options.cells) == 'function') {
                $.each(cells, function (k, v) {
                    var coords = $(v).prop('id').split('-');
                    options.cells($(v), coords[0], coords[1]);
                });
            } else if (typeof(options.table) == 'function') {
                $.each(cells, function (k, v) {
                    var coords = $(v).prop('id').split('-');
                    options.table($(main), $(v), coords[0], coords[1], $.fn.jexcel.defaults[id].data[coords[1]][coords[0]], $.fn.jexcel('getColumnNameFromId', [coords[0], coords[1]]));
                });
            }
        }
    },

    /**
     * Open the editor
     * 
     * @param object cell
     * @return void
     */
    openEditor : function(cell, empty, e) {
        // Id
        var id = $(this).prop('id');

        // Main
        var main = $(this);

        // Options
        var options = $.fn.jexcel.defaults[id];

        // Get cell position
        var position = $(cell).prop('id').split('-');

        // Readonly
        if ($(cell).hasClass('readonly') == true) {
            // Do nothing
        } else {
            // Holder
            $.fn.jexcel.edition = $(cell).html();
            $.fn.jexcel.editionValue = $(this).jexcel('getValue', $(cell));

            // If there is a custom editor for it
            if (options.columns[position[0]].editor) {
                // Keep the current value
                $(cell).addClass('edition');

                // Custom editors
                options.columns[position[0]].editor.openEditor(cell, empty, e);
            } else {
                // Native functions
                if (options.columns[position[0]].type == 'checkbox' || options.columns[position[0]].type == 'radio' || options.columns[position[0]].type == 'hidden') {
                    // Get value
                    var value = $(this).jexcel('getValue', $(cell)) == 1 ? false : true;
                    // Update value
                    $(this).jexcel('setValue', $(cell), value);
                } else if (options.columns[position[0]].type == 'dropdown' || options.columns[position[0]].type == 'autocomplete') {
                    // Keep the current value
                    $(cell).addClass('edition');

                    // Get current value
                    var value = $(cell).find('input').val();

                    // Create dropdown
                    if (typeof(options.columns[position[0]].filter) == 'function') {
                        var source = options.columns[position[0]].filter($(this), $(cell), position[0], position[1], options.columns[position[0]].source);
                    } else {
                        var source = options.columns[position[0]].source;
                    }

                    if ($.fn.jdropdown) {
                        // Create dropdown
                        var editor = document.createElement('div');
                        $(editor).jdropdown({
                            data: source,
                            multiple: options.columns[position[0]].multiple ? true : false,
                            autocomplete: options.columns[position[0]].autocomplete || options.columns[position[0]].type == 'autocomplete' ? true : false,
                            opened:true,
                            value: options.columns[position[0]].multiple ? value.split(';') : value,
                            onclose:function() {
                                $(main).jexcel('closeEditor', $(cell), true);
                            }
                        });

                        // Editor configuration
                        $(editor).css('width', $(cell).width());
                        $(editor).css('height', $(cell).height());
                    } else {
                        if (options.columns[position[0]].type == 'autocomplete') {
                            // This is deprected. You should include jdropdown to render the autocomplete properly 
                            console.error('This is deprected. You should include jdropdown to render the autocomplete properly');
                        } else {
                            // Create HTML content
                            var editor = document.createElement('select');
                            for (i = 0; i < source.length; i++) {
                                if (typeof(source[i]) == 'object') {
                                    k = source[i].id;
                                    v = source[i].name;
                                } else {
                                    k = source[i];
                                    v = source[i];
                                }
                                var option = document.createElement('option');
                                $(option).val(k);
                                $(option).html(v);
                                $(editor).append(option);
                            }
                            // Editor configuration
                            $(editor).css('width', $(cell).width());
                            $(editor).css('height', $(cell).height());
                            // Blur
                            $(editor).blur(function () {
                                $(main).jexcel('closeEditor', $(this).parent(), true);
                            });
                            // On change
                            $(editor).change(function () {
                                $(main).jexcel('closeEditor', $(this).parent(), true);
                            });
                            // Focus
                            $(editor).focus();
                            // Set value
                            if (value) {
                                $(editor).val(value);
                            }
                        }
                    }

                    // Append editor to the cell
                    $(cell).html(editor);
                } else if (options.columns[position[0]].type == 'calendar') {
                    $(cell).addClass('edition');
                    // Get current value
                    var value = $(cell).find('input').val();
                    // Basic editor
                    var editor = document.createElement('input');
                    $(editor).prop('class', 'editor');
                    $(editor).css('width', $(cell).width());
                    $(editor).css('height', $(cell).height());
                    $(editor).val($(cell).text());
                    $(cell).html(editor);
                    $(editor).focus();

                    options.columns[position[0]].options.onclose = function () {
                        $(main).jexcel('closeEditor', $(cell), true);
                    }

                    // Current value
                    $(editor).jcalendar(options.columns[position[0]].options);
                    $(editor).jcalendar('open', value);
                } else {
                    // Get cell proprieties
                    var editorWidth = $(cell).width();
                    var editorHeight = $(cell).innerHeight();

                    // Keep the current value
                    $(cell).addClass('edition');

                    if (options.wordWrap == true || options.columns[position[0]].wordWrap == true) {
                        var input = $(cell).find('textarea');
                    } else {
                        var input = $(cell).find('input');
                    }

                    // Get content
                    if ($(input).length) {
                        var html = $(input).val();
                    } else {
                        var html = $(cell).html();
                    }

                    // Basic editor
                    if (options.wordWrap == true || options.columns[position[0]].wordWrap == true) {
                        var editor = document.createElement('textarea');
                    } else {
                        var editor = document.createElement('input');
                    }

                    $(editor).prop('class', 'editor');
                    $(editor).css('width', editorWidth);
                    $(editor).css('min-height', editorHeight);
                    $(cell).html(editor);

                    // Bind mask
                    if (options.columns[position[0]].mask) {
                        if (! $.fn.masked) {
                            console.error('Jexcel: it was not possible to load the mask plugin.');
                        } else {
                            $(editor).mask(options.columns[position[0]].mask, options.columns[position[0]].options)
                        }
                    }

                    // Current value
                    $(editor).focus();
                    if (! empty) {
                        $(editor).val(html);
                    }

                    // Close editor handler
                    $(editor).blur(function () {
                        $(main).jexcel('closeEditor', $(this).parent(), true);
                    });
                }
            }
        }
    },

    /**
     * Close the editor and save the information
     * 
     * @param object cell
     * @param boolean save
     * @return void
     */
    closeEditor : function(cell, save) {
        // Remove edition mode mark
        $(cell).removeClass('edition');

        // Id
        var id = $(this).prop('id');

        // Options
        var options = $.fn.jexcel.defaults[id];

        // Cell identification
        var position = $(cell).prop('id').split('-');

        // Get cell properties
        if (save == true) {
            // If custom editor
            if (options.columns[position[0]].editor) {
                // Custom editor
                value = options.columns[position[0]].editor.closeEditor(cell, save);
            } else {
                // Native functions
                if (options.columns[position[0]].type == 'checkbox' || options.columns[position[0]].type == 'radio' || options.columns[position[0]].type == 'hidden') {
                    // Do nothing
                } else if (options.columns[position[0]].type == 'dropdown' || options.columns[position[0]].type == 'autocomplete') {
                    // Get value
                    if ($.fn.jdropdown) {
                        var value = $(cell).find('.jdropdown').jdropdown('getValue');
                    } else {
                        if (options.columns[position[0]].type == 'autocomplete') {
                            console.error('This is deprected. You should include jdropdown to render the autocomplete properly');
                        } else {
                            var value = $(cell).find('select').val();
                        }
                    }
                } else if (options.columns[position[0]].type == 'calendar') {
                    var value = $(cell).find('.jcalendar_value').val();
                } else if (options.columns[position[0]].type == 'numeric') {
                    var value = $(cell).find('.editor').val();
                    if (value.substr(0,1) != '=') {
                        if (value == '' && options.columns[position[0]].allowEmpty) {
                            value = '';
                        } else {
                            value = Number(value) || 0;
                        }
                    }
                } else {
                    // Get content
                    var value = $(cell).find('.editor').val();
                }
            }

            // Update values
            var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

            // Ignore changes if the value is the same
            if ($.fn.jexcel.defaults[id].data[position[1]][position[0]] == value) {
                // Disabled events and history
                $.fn.jexcel.ignoreHistory = true;
            }

            // Get value from column and set the default
            $.fn.jexcel.defaults[id].data[position[1]][position[0]] = value;

            // Update cell
            $(this).jexcel('updateCells', [{
                col: position[0],
                row: position[1],
                cell: $(cell),
                newValue: value,
                oldValue: $.fn.jexcel.editionValue
            }]);

            // Restore events and history flag
            $.fn.jexcel.ignoreHistory = ignoreHistory;
        } else {
            if (options.columns[position[0]].editor) {
                // Custom editor
                options.columns[position[0]].editor.closeEditor(cell, save);
            } else if (options.columns[position[0]].type == 'calendar') {
                // Do nothing - calendar will be closed without keeping the current value
            } else {
                if (options.columns[position[0]].type == 'dropdown' || options.columns[position[0]].type == 'autocomplete') {
                    $(cell).find('.jdropdown').jdropdown('close');
                }

                // Restore value
                $(cell).html($.fn.jexcel.edition);

                // Finish temporary edition
                $.fn.jexcel.edition = null;
            }
        }
    },

    /**
     * Get the cell object
     * 
     * @param object cell
     * @return string value
     */
    getCell : function(cell) {
        // Convert in case name is excel liked ex. A10, BB92
        cell = $(this).jexcel('getIdFromColumnName', cell);
        // Get object based on a string ex. 12-1, 13-3
        cell = $(this).find('[id=' + cell +']');

        return cell;
    },

    /**
     * Get the value from a cell
     * 
     * @param object cell
     * @return string value
     */
    getValue : function(cell) {
        var value = null;

        // If is a string get the cell object
        if (typeof(cell) != 'object') {
            // Get cell excel like A11, B99, etc
            cell = $(this).jexcel('getCell', cell);
        }

        // If column exists
        if ($(cell).length) {
            // Id
            var id = $(this).prop('id');

            // Global options
            var options = $.fn.jexcel.defaults[id];

            // Configuration
            var position = $(cell).prop('id').split('-');

            // Get value based on the type
            if (options.columns[position[0]].editor) {
                // Custom editor
                value = options.columns[position[0]].editor.getValue(cell);
            } else {
                // Native functions
                if (options.columns[position[0]].type == 'checkbox' || options.columns[position[0]].type == 'radio') {
                    // Get checkbox value
                    value = $(cell).find('input').is(':checked') == true ? '1' : '0';
                } else if (options.columns[position[0]].type == 'dropdown' || options.columns[position[0]].type == 'autocomplete' || options.columns[position[0]].type == 'calendar') {
                    // Get value
                    value = $(cell).find('input').val();
                } else if (options.columns[position[0]].type == 'currency') {
                    value = $(cell).html().replace( /\D/g, '');
                } else {
                    // Get default value
                    value = $(cell).find('input');
                    if ($(value).length) {
                        value = $(value).val(); 
                    } else {
                        value = $(cell).html();
                    }
                }
            }
        }

        return value;
    },

    /**
     * Set a cell value
     * 
     * @param object cell destination cell
     * @param object value value
     * @return void
     */
    setValue : function(cell, value, force) {
        // Check request type, object or column name
        if (typeof(cell) != 'object') {
            // Get cell identification
            var cellName = $(this).jexcel('getIdFromColumnName', cell);
            // Get cell object based on excel names A1, A2, etc...
            cell = $(this).jexcel('getCell', cell);
        }

        // If column exists
        if ($(cell).length) {
            // Id
            var id = $(this).prop('id');

            // Main object
            var main = $(this);

            // Records to be updated
            var records = [];

            // Go throw all cells
            $.each(cell, function(k, v) {
                // Cell position
                var position = $(v).prop('id').split('-');

                // Update cell
                records.push({
                    col: position[0],
                    row: position[1],
                    cell: $(v),
                    newValue: value,
                    oldValue: $.fn.jexcel.defaults[id].data[position[1]][position[0]],
                });
            });

            // Force change for readonly columns
            force = force ? true : false;

            // Update cells
            $(this).jexcel('updateCells', records, force);

            return true;
        } else {
            return false;
        }
    },

    /**
     * Get the label from a cell
     * 
     * @param object cell
     * @return string value
     */
    getText : function(cell) {
        var value = null;

        // If is a string get the cell object
        if (typeof(cell) != 'object') {
            // Get cell excel like A11, B99, etc
            cell = $(this).jexcel('getCell', cell);
        }

        // If column exists
        if ($(cell).length) {
            // Id
            var id = $(this).prop('id');

            // Global options
            var options = $.fn.jexcel.defaults[id];

            // Configuration
            var position = $(cell).prop('id').split('-');

            // Get value based on the type
            if (options.columns[position[0]].editor) {
                // Custom editor
                value = options.columns[position[0]].editor.getText(cell);
            } else {
                // Native functions
                if (options.columns[position[0]].type == 'checkbox' || options.columns[position[0]].type == 'radio') {
                    // Get checkbox value
                    value = $(cell).find('input').is(':checked') == true ? '1' : '0';
                } else {
                    value = $(cell).text();
                }
            }
        }

        return value;
    },

    /**
     * Update cells with no history and events
     * 
     * @param object destination cells
     * @param bool ignoreHistory - keep cell change out of the undo/redo history
     * @return void
     */
    loadCells : function(cells, force) {
        // Update values
        var ignoreEvents = $.fn.jexcel.ignoreEvents ? true : false;
        var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

        // Disabled events and history
        $.fn.jexcel.ignoreEvents = true;
        $.fn.jexcel.ignoreHistory = true;

        // Update all records
        $(this).jexcel('updateCells', cells, force);

        // Restore events and history flag
        $.fn.jexcel.ignoreEvents = ignoreEvents;
        $.fn.jexcel.ignoreHistory = ignoreHistory;
    },

    /**
     * Update cells content
     * 
     * @param object destination cells
     * @param bool ignoreHistory - keep cell change out of the undo/redo history
     * @return void
     */
    updateCells : function(cells, force) {
        // Id
        var id = $(this).prop('id');

        // Main object
        var main = $(this);

        // Global options
        var options = $.fn.jexcel.defaults[id];

        // Update cells
        $.each(cells, function (k, v) {
            // Get cell
            v.cell = $(main).find('#' + v.col + '-' + v.row); // TODO: check if v.cell is necessary

            // Before Change
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof(options.onbeforechange) == 'function') {
                    options.onbeforechange(main, $(v.cell), v.oldValue, v.newValue);
                }
            }

            // Force
            force = force ? true : false;

            // Update
            $(main).jexcel('updateCell', v, force);

            // Change
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof(options.onchange) == 'function') {
                    options.onchange(main, $(v.cell), v.newValue, v.oldValue);
                }
            }
        });

        // After changes
        $(this).jexcel('afterChange');

        // Keeping history of changes
        $(this).jexcel('setHistory', cells);
    },

    /**
     * Update cell content
     * 
     * @param object cell
     * @return void
     */
    updateCell : function(v, force) {
        // Id
        var id = $(this).prop('id');

        // Global options
        var options = $.fn.jexcel.defaults[id];

        // Cell identification
        var position = [v.col, v.row];

        // Remove formula error
        $(v.cell).removeClass('error');

        // Value
        value = '' + v.newValue;

        // Changing value depending on the column type
        if ($(v.cell).hasClass('readonly') == true && force == false) {
            // Do nothing
        } else {
            if (options.columns[position[0]].editor) {
                // Custom editor
                options.columns[position[0]].editor.setValue($(v.cell), value);
            } else {
                // Native functions
                if (options.columns[position[0]].type == 'checkbox' || options.columns[position[0]].type == 'radio') {
                    if (value == 1 || value == true || value == 'true') {
                        $(v.cell).find('input').prop('checked', true);
                    } else {
                        $(v.cell).find('input').prop('checked', false);
                    }
                } else if (options.columns[position[0]].type == 'dropdown' || options.columns[position[0]].type == 'autocomplete') {
                    // Dropdown and autocompletes
                    key = '';
                    val = '';
                    if (value) {
                        // Realtime combo to get correct labels
                        var combo = [];
                        var source = options.columns[position[0]].source;

                        for (num = 0; num < source.length; num++) {
                            if (typeof(source[num]) == 'object') {
                                combo[source[num].id] = source[num].name;
                            } else {
                                combo[source[num]] = source[num];
                            }
                        }

                        // Garante single multiple compatibily
                        var values = value.split(';')

                        for (var values_i = 0; values_i < values.length; values_i++) {
                            if (combo[values[values_i]]) {
                                // Keys
                                if (key) {
                                    key += ';'
                                }
                                key += values[values_i];
                                // Value
                                if (val) {
                                    val += '; '
                                }
                                val += combo[values[values_i]];
                            }
                        }
                    }

                    if (! val) {
                        val = '&nbsp';
                    }

                    $(v.cell).html('<input type="hidden" value="' +  key + '">' + val + '<span class="jexcel_arrow"><span id="jexcel_arrow"></span></span>');
                } else if (options.columns[position[0]].type == 'calendar') {
                    val = '';
                    if (value != 'undefined') {
                        // Formated?
                        date = $.fn.jcalendar('fromFormatted', value, options.columns[position[0]].options.format);
                        if (date) {
                            value = date;
                        }

                        val = $.fn.jcalendar('label', value, options.columns[position[0]].options.format);
                    } else {
                        val = '';
                    }

                    if (! val) {
                        value = '';
                    }

                    $(v.cell).html('<input type="hidden" value="' + value + '">' + val);
                } else {
                    val = value ? value : '';

                    // Is this value a formula
                    if (value.substr(0,1) == '=') {
                        $(this).jexcel('updateFormulas', v.col + '-' + v.row, value);

                        $(v.cell).html('<input type="hidden" value="' + value + '">' + val);
                    } else {
                        // This is not a formula
                        if (options.columns[position[0]].type == 'numeric') {
                            $(v.cell).html('<input type="hidden" value="' + value + '">' + val);
                        } else {
                            $(v.cell).html(val);
                        }
                    }
                }
            }

            // Get value from column and set the default
            $.fn.jexcel.defaults[id].data[position[1]][position[0]] = value;

            // Global variables
            var i = $.fn.jexcel('getColumnName', v.col)  + (parseInt(v.row) + 1);

            if (value == Number(value)) {
                $.fn.jexcel.defaults[id].values[i] = Number(value);
            } else {
                $.fn.jexcel.defaults[id].values[i] = value;
            }

            // Local spreasheet values
            window[i] = new $.fn.jexcel.factory(i, $.fn.jexcel.defaults[id].values[i]);

            // Events?
            if ($.fn.jexcel.ignoreEvents != true) {
                // Update own cell
                if (value.substr(0,1) == '=') {
                    $(this).jexcel('executeFormula', v.col + '-' + v.row);
                }

                // Update related cells
                if (Object.keys($.fn.jexcel.defaults[id].formula).length) {
                    $(this).jexcel('formula', v.col + '-' + v.row);
                }
            }
        }
    },

    /**
     * Update the cells selection
     * 
     * @param object o cell origin
     * @param object d cell destination
     * @return void
     */
    updateSelection : function(o, d, origin) {
        // Main table
        var main = $(this);

        // Id
        var id = $(this).prop('id');

        // Cells
        var rows = $(this).find('tbody tr');
        var cells = $(this).find('tbody td');
        var header = $(this).find('.jexcel_headers td');

        // Keep previous selection status
        var previousStatus = ($(this).find('.highlight').length > 0) ? true : false;
        var currentStatus = false;

        // Remove highlight
        $(cells).removeClass('highlight');
        $(cells).removeClass('highlight-left');
        $(cells).removeClass('highlight-right');
        $(cells).removeClass('highlight-top');
        $(cells).removeClass('highlight-bottom');

        // Update selected column
        $(header).removeClass('selected');
        $(rows).removeClass('selected');

        // Origin & Destination
        if (o && d) {
            // Define coordinates
            or = $(o).prop('id').split('-');
            de = $(d).prop('id').split('-');

            if (parseInt(or[0]) < parseInt(de[0])) {
                var px = parseInt(or[0]);
                var ux = parseInt(de[0]);
            } else {
                var px = parseInt(de[0]);
                var ux = parseInt(or[0]);
            }

            if (parseInt(or[1]) < parseInt(de[1])) {
                var py = parseInt(or[1]);
                var uy = parseInt(de[1]);
            } else {
                var py = parseInt(de[1]);
                var uy = parseInt(or[1]);
            }

            // Limits
            var borderLeft = null;
            var borderRight = null;
            var borderTop = py;
            var borderBottom = uy;

            // Redefining styles
            for (var i = px; i <= ux; i++) {
                for (var j = py; j <= uy; j++) {
                    // Select markers
                    $(this).find('#' + i + '-' + j).addClass('highlight');

                    // Row and column headers
                    $(this).find('#col-' + i).addClass('selected');
                    $(this).find('#row-' + j).addClass('selected');
                }

                // Right limits
                if ($.fn.jexcel.defaults[id].columns[i].type != 'hidden') {
                    if (borderLeft == null) {
                        borderLeft = i;
                    }
                    borderRight = i;
                }
            }

            // Create borders
            for (var i = borderLeft; i <= borderRight; i++) {
                // Top border
                $(this).find('#' + i + '-' + borderTop).addClass('highlight-top');
                // Bottom border
                $(this).find('#' + i + '-' + borderBottom).addClass('highlight-bottom');
            }

            for (var j = borderTop; j <= borderBottom; j++) {
                // Left border
                $(this).find('#' + borderLeft + '-' + j).addClass('highlight-left');
                // Right border
                $(this).find('#' + borderRight + '-' + j).addClass('highlight-right');
            }

            // Get current selection status
            var currentStatus = true;
        }

        if ($.fn.jexcel.ignoreEvents != true) {
            if ($.fn.jexcel.defaults[id].onblur) {
                if (typeof($.fn.jexcel.defaults[id].onblur) == 'function') {
                    if (previousStatus == true && currentStatus == false) {
                        $.fn.jexcel.defaults[id].onblur($(this));
                    }
                }
            }

            if ($.fn.jexcel.defaults[id].onfocus) {
                if (typeof($.fn.jexcel.defaults[id].onfocus) == 'function') {
                    if (previousStatus == false && currentStatus == true) {
                        $.fn.jexcel.defaults[id].onfocus($(this));
                    }
                }
            }

            if (currentStatus == true) {
                // Events
                if ($.fn.jexcel.defaults[id].onselection) {
                    if (typeof($.fn.jexcel.defaults[id].onselection) == 'function') {
                        $.fn.jexcel.defaults[id].onselection($(this), o, d, origin);
                    }
                }
            }
        }

        // Keep real selection : different from visual selection
        if (o && d) {
            $.fn.jexcel.selection = [o, d];
        } else {
            $.fn.jexcel.selection = null;
        }

        // Find corner cell
        setTimeout(function() {
            $(main).jexcel('updateCornerPosition');
        }, 0);
    },

    /**
     * Reset all selection from the current table
     */
    resetSelection : function() {
        // Remove selection from any other jexcel if applicable
        if ($.fn.jexcel.current) {
            $('#' + $.fn.jexcel.current).jexcel('updateSelection');
        }

        // Hide corner
        $('.jexcel_corner').css('left', '-200px');

        // Reset controls
        $.fn.jexcel.selectedRow = null;
        $.fn.jexcel.selectedCell = null;
        $.fn.jexcel.selectedHeader = null;
    },

    /**
     * Get the selection parameters
     * 
     * @return array [first cell, last cell]
     */
    getSelection : function() {
        // Return selection
        if ($.fn.jexcel.current) {
            return $.fn.jexcel.selection;
        }
    },

    /**
     * Current the selected cursor cell
     * 
     * @return object cell
     */
    getCellCursor : function() {
        if ($.fn.jexcel.current) {
            return $.fn.jexcel.selectedCell;
        }
    },

    /**
     * Get all the selected cells in the current table
     * 
     * @return object cells
     */ 
    getSelectedCells : function() {
        // Get selected cells
        if ($.fn.jexcel.current) {
            return $('#' + $.fn.jexcel.current).find('tbody td.highlight');
        }
    },

    /**
     * Update the cells move data
     * 
     * @param object o cell current
     * @return void
     */
    updateCornerSelection : function(current) {
        // Main table
        var main = $(this);

        // Remove selection
        var cells = $(this).find('tbody td');
        $(cells).removeClass('selection');
        $(cells).removeClass('selection-left');
        $(cells).removeClass('selection-right');
        $(cells).removeClass('selection-top');
        $(cells).removeClass('selection-bottom');

        // Get selection
        var selection = $(this).find('tbody td.highlight');

        // Get elements first and last
        var s = $(selection[0]).prop('id').split('-');
        var d = $(selection[selection.length - 1]).prop('id').split('-');

        // Get current
        var c = $(current).prop('id').split('-');

        // Vertical copy
        if (c[1] > d[1] || c[1] < s[1]) {
            // Vertical
            var px = parseInt(s[0]);
            var ux = parseInt(d[0]);
            if (parseInt(c[1]) > parseInt(d[1])) {
                var py = parseInt(d[1]) + 1;
                var uy = parseInt(c[1]);
            } else {
                var py = parseInt(c[1]);
                var uy = parseInt(s[1]) - 1;
            }
        } else if (c[0] > d[0] || c[0] < s[0]) {
            // Horizontal copy
            var py = parseInt(s[1]);
            var uy = parseInt(d[1]);
            if (parseInt(c[0]) > parseInt(d[0])) {
                var px = parseInt(d[0]) + 1;
                var ux = parseInt(c[0]);
            } else {
                var px = parseInt(c[0]);
                var ux = parseInt(s[0]) - 1;
            }
        }

        for (j = py; j <= uy; j++) {
            for (i = px; i <= ux; i++) {
                $(this).find('#' + i + '-' + j).addClass('selection');
                $(this).find('#' + i + '-' + py).addClass('selection-top');
                $(this).find('#' + i + '-' + uy).addClass('selection-bottom');
                $(this).find('#' + px + '-' + j).addClass('selection-left');
                $(this).find('#' + ux + '-' + j).addClass('selection-right');
            }
        }
    },

    /**
     * Update corner position
     * 
     * @return void
     */
    updateCornerPosition : function() {
        // Highlighted cells
        var cells = $(this).find('.highlight');

        // If any selected cells
        if ($(cells).length) { 
            // Get last cell
            var corner = $(cells).last();

            // Get the position of the corner helper
            var t = parseInt($(corner).offset().top) + $(corner).parent().innerHeight() - 3;
            var l = parseInt($(corner).offset().left) + $(corner).innerWidth() - 2;

            // Place the corner in the correct place
            $('.jexcel_corner').css('top', t);
            $('.jexcel_corner').css('left', l);

            // Hide the corner in case is out of the range
            var docViewTop = $(this).offset().top;
            var docViewBottom = docViewTop + $(this).height();
            var elemTop = t;
            var elemBottom = t;

            if (! ((elemBottom <= docViewBottom) && (elemTop >= docViewTop))) {
                $('.jexcel_corner').css('top', -200);
                $('.jexcel_corner').css('left', -200);
            }

            var docViewLeft = $(this).offset().left;
            var docViewRight = docViewLeft + $(this).width();
            var elemLeft = l;
            var elemRight = l;

            if (! ((elemRight <= docViewRight - 5) && (elemLeft >= docViewLeft))) {
                $('.jexcel_corner').css('top', -200);
                $('.jexcel_corner').css('left', -200);
            }
        }
    },

    /**
     * Get the whole table data
     * 
     * @param integer row number
     * @return string value
     */
    getData : function(highlighted) {
        // Control vars
        var dataset = [];
        var px = 0;
        var py = 0;

        var id = $(this).prop('id');

        // Column and row length
        var x = $.fn.jexcel.defaults[id].columns.length;
        var y = $.fn.jexcel.defaults[id].data.length;

        // Go through the columns to get the data
        for (j = 0; j < y; j++) {
            px = 0;
            for (i = 0; i < x; i++) {
                // Cell
                cell = $(this).find('#' + i + '-' + j);

                // Cell selected or fullset
                if (! highlighted || $(cell).hasClass('highlight')) {
                    // Get value
                    if (! dataset[py]) {
                        dataset[py] = [];
                    }
                    dataset[py][px] = $(this).jexcel('getValue', $(cell));
                    px++;
                }
            }
            if (px > 0) {
                py++;
            }
        }

       return dataset;
    },

    /**
     * Get a row data by row number
     * 
     * @param integer rowNumber
     * @return array
     */
    getRowData : function(rowNumber) {
        // Main instance
        var main = $(this);
        // Data to be returned
        var dataset = [];

        // Row number
        if (! parseInt(rowNumber)) {
            rowNumber = 0;
        } else if (rowNumber == -1) {
            // Id
            var id = $(this).prop('id');
            // Get the main object configuration
            var options = $.fn.jexcel.defaults[id];
            // Last row
            rowNumber = options.data.length - 1;
        }

        // Selected row
        var cells = $(this).find('.jexcel-content .r' + rowNumber);

        // Go through the columns to get the data
        $.each(cells, function(k, v) {
            dataset.push($(main).jexcel('getValue', $(v)));
        });

        return dataset;
    },

    /**
     * Get a column data by column number
     * 
     * @param integer columnNumber
     * @return array
     */
    getColumnData : function(columnNumber) {
        // Main instance
        var main = $(this);
        // Data to be returned
        var dataset = [];

        // Row number
        if (! parseInt(columnNumber)) {
            columnNumber = 0;
        } else if (columnNumber == -1) {
            // Id
            var id = $(this).prop('id');
            // Get the main object configuration
            var options = $.fn.jexcel.defaults[id];
            // Last column
            columnNumber = options.columns.length - 1;
        }

        // Selected row
        var cells = $(this).find('.jexcel-content .c' + columnNumber);

        // Go through the columns to get the data
        $.each(cells, function(k, v) {
            dataset.push($(main).jexcel('getValue', $(v)));
        });

        return dataset;
    },

    /**
     * Copy method
     * 
     * @param bool highlighted - Get only highlighted cells
     * @param delimiter - \t default to keep compatibility with excel
     * @return string value
     */
    copy : function(highlighted, delimiter, returnData) {
        // Id
        var id = $(this).prop('id');

        if (! delimiter) {
            delimiter = "\t";
        }

        // Column and row length
        var col = [];
        var row = [];
        var x = $.fn.jexcel.defaults[id].data[0].length
        var y = $.fn.jexcel.defaults[id].data.length
        var tmp = '';

        // Reset container
        var style = [];

        // Go through the columns to get the data
        for (var j = 0; j < y; j++) {
            col = [];

            for (var i = 0; i < x; i++) {
                // Get cell
                cell = $(this).find('#' + i + '-' + j);

                // If cell is highlighted
                if (! highlighted || $(cell).hasClass('highlight')) {
                    // Get value
                    val = $(this).jexcel('getValue', $(cell));
                    if (val.match && (val.match(/,/g) || val.match(/\n/) || val.match(/\"/))) {
                        // Scape double quotes
                        val = val.replace(new RegExp('"', 'g'), '""');
                        val = '"' + val + '"'; 
                    }
                    col.push(val);
                    // Get style
                    tmp = $(cell).attr('style');
                    style.push(tmp ? tmp : '');
                }
            }

            if (col.length) {
                row.push(col.join(delimiter));
            }
        }

        // Final string
        var str = row.join("\n");

        // Create a hidden textarea to copy the values
        if (! returnData) {
            txt = $('.jexcel_textarea');
            $(txt).val(str);
            $(txt).select();
            document.execCommand("copy");
        }

        // Keep non visible information
        $.fn.jexcel.hash = $(this).jexcel('hash', str); 
        $.fn.jexcel.style = style;

        return str;
    },

    /**
     * jExcel cut method
     */ 
    cut : function () {
        // Reset cells
        var cells = $(this).find('.highlight');

        // If cell exists
        if ($(cells).length > 0) {
            var records = [];
            // Save current cell to the history
            $.each(cells, function(k, v) {
                // Get cell information
                var position = $(v).prop('id').split('-');

                // Keep history
                records.push({
                    col: position[0],
                    row: position[1],
                    cell: $(v),
                    newValue: '',
                    oldValue: $.fn.jexcel.defaults[$.fn.jexcel.current].data[position[1]][position[0]],
                });
            });

            // Copy data
            $(this).jexcel('copy', true);

            // Reset cells
            $(this).jexcel('updateCells', records);
        }
    },

    /**
     * jExcel paste method
     * 
     * @param integer row number
     * @return string value
     */
    paste : function(cell, data) {
        // Id
        var id = $(this).prop('id');

        // Paste filter
        if (typeof($.fn.jexcel.defaults[id].onbeforepaste) == 'function') {
            data = $.fn.jexcel.defaults[id].onbeforepaste(data);
        }

        // Possible style
        var hash = $(this).jexcel('hash', data);
        var style = (hash == $.fn.jexcel.hash) ? $.fn.jexcel.style : null;

        // Parse paste
        data = $(this).jexcel('parseCSV', data, "\t");

        // Initial position
        var position = $(cell).prop('id');

        if (position) {
            position = position.split('-');
            var x = parseInt(position[0]);
            var y = parseInt(position[1]);

            // Automatic adding new rows when the copied data is larger then the table
            if (y + data.length > $.fn.jexcel.defaults[id].data.length) {
                $(this).jexcel('insertRow', y + data.length - $.fn.jexcel.defaults[id].data.length);
            }
            // Automatic adding new columns when the copied data is larger then the table
            if (data[0]) {
                row = data[0];
                if (x + row.length > $.fn.jexcel.defaults[id].data[y].length) {
                    $(this).jexcel('insertColumn', x + row.length - $.fn.jexcel.defaults[id].data[y].length);
                }
            }

            // Records
            var records = []; 

            // Index
            var styleIndex = 0;

            // Go through the columns to get the data
            for (j = 0; j < data.length; j++) {
                // Explode column values
                row = data[j];
                for (i = 0; i < row.length; i++) {
                    // Get cell
                    cell = $(this).find('#' + (parseInt(i) + parseInt(x))  + '-' + (parseInt(j) + parseInt(y)));

                    // If cell exists
                    if ($(cell).length > 0) {
                        // Keep cells history
                        records.push({
                            col: (parseInt(i) + parseInt(x)),
                            row: (parseInt(j) + parseInt(y)),
                            cell: $(cell),
                            newValue: row[i],
                            oldValue: $.fn.jexcel.defaults[id].data[(parseInt(j) + parseInt(y))][(parseInt(i) + parseInt(x))],
                        });
                    }

                    // Style
                    if (style) {
                        $(cell).attr('style', style[styleIndex]);
                    }

                    styleIndex++;
                }
            }

            // Select the new cells
            $(this).jexcel('updateSelection', records[0].cell, records[records.length-1].cell);

            // Save history
            if (records.length > 0) {
                // Update new values
                $(this).jexcel('updateCells', records);

                // Paste event
                if (typeof($.fn.jexcel.defaults[id].onpaste) == 'function') {
                    data = $.fn.jexcel.defaults[id].onpaste(records);
                }
            }
        }
    },

    /**
     * Based on script by Ben Nadel
     */
    parseCSV : function(str, delimiter) {
        // Remove last line break
        str = str.replace(/\r?\n$|\r$|\n$/g, "");
        // Last caracter is the delimiter
        if (str.charCodeAt(str.length-1) == 9) {
            str += "\0";
        }
        // user-supplied delimeter or default comma
        delimiter = (delimiter || ",");

        var arr = [];
        var quote = false;  // true means we're inside a quoted field
        // iterate over each character, keep track of current row and column (of the returned array)
        for (var row = 0, col = 0, c = 0; c < str.length; c++) {
            var cc = str[c], nc = str[c+1];
            arr[row] = arr[row] || [];
            arr[row][col] = arr[row][col] || '';

            // If the current character is a quotation mark, and we're inside a quoted field, and the next character is also a quotation mark, add a quotation mark to the current column and skip the next character
            if (cc == '"' && quote && nc == '"') { arr[row][col] += cc; ++c; continue; }  

            // If it's just one quotation mark, begin/end quoted field
            if (cc == '"') { quote = !quote; continue; }

            // If it's a comma and we're not in a quoted field, move on to the next column
            if (cc == delimiter && !quote) { ++col; continue; }

            // If it's a newline (CRLF) and we're not in a quoted field, skip the next character and move on to the next row and move to column 0 of that new row
            if (cc == '\r' && nc == '\n' && !quote) { ++row; col = 0; ++c; continue; }

            // If it's a newline (LF or CR) and we're not in a quoted field, move on to the next row and move to column 0 of that new row
            if (cc == '\n' && !quote) { ++row; col = 0; continue; }
            if (cc == '\r' && !quote) { ++row; col = 0; continue; }

            // Otherwise, append the current character to the current column
            arr[row][col] += cc;
        }
        return arr;
    },

    /**
     * Insert a new column
     * 
     * @param  mixed - num of columns to be added or data to be added in one single column
     * @param  object properties - column properties
     * @param  int numColumns - number of columns to be created
     * @return void
     */
    insertColumn : function (mixed, properties, columnNumber) {
        // Id
        var id = $(this).prop('id');

        // Get the main object configuration
        var options = $.fn.jexcel.defaults[id];

        // Configuration
        if (options.allowInsertColumn == true) {
            // Records
            var records = [];

            // Data to be insert
            var data = [];

            // The insert could be lead by number of rows or the array of data
            if (mixed > 0) {
                numOfColumns = mixed;
            } else {
                numOfColumns = 1;

                if (mixed) {
                    data = mixed;
                }
            }

            // Current column number
            var lastColumn = options.colHeaders.length;

            // Append direction
            var direction = 0;

            // Confirm position
            if (columnNumber == undefined || columnNumber >= parseInt(lastColumn)) {
                // Default is the end of the table
                columnNumber = parseInt(lastColumn);

                // Direction
                var direction = 1;
            }

            // Create default properties
            if (! properties) {
                properties = [];
            }
            if (! properties.columns) {
                properties.columns = [];
            }
            if (! properties.colWidths) {
                properties.colWidths = [];
            }
            if (! properties.colHeaders) {
                properties.colHeaders = [];
            }
            if (! properties.colAlignments) {
                properties.colAlignments = [];
            }
            if (! properties.colHeaderClasses) {
                properties.colHeaderClasses = [];
            }

            for (var i = 0; i < numOfColumns; i++) {
                if (! properties.columns[i]) {
                    properties.columns[i] = { type:'text' };
                }
                if (! properties.colWidths[i]) {
                    properties.colWidths[i] = 50;
                }
                if (! properties.colHeaders[i]) {
                    properties.colHeaders[i] = '';
                }
                if (! properties.colAlignments[i]) {
                    properties.colAlignments[i] = 'center';
                }
                if (! properties.colHeaderClasses[i]) {
                    properties.colHeaderClasses[i] = '';
                }
            }

            // Keep current configuration
            var currentColumns = options.columns.splice(columnNumber);
            var currentColWidths = options.colWidths.splice(columnNumber);
            var currentColHeaders = options.colHeaders.splice(columnNumber);
            var currentColAlignments = options.colAlignments.splice(columnNumber);
            var currentColHeaderClasses = options.colHeaderClasses.splice(columnNumber);

            // Adding headers configuration
            i = 0;

            for (var col = columnNumber; col < (numOfColumns + columnNumber); col++) {
                // Adding the column properties to the main property holder
                options.columns[col] = properties.columns[i];
                options.colWidths[col] = properties.colWidths[i];
                options.colHeaders[col] = properties.colHeaders[i];
                options.colAlignments[col] = properties.colAlignments[i];
                options.colHeaderClasses[col] = properties.colHeaderClasses[i];

                if (! options.columns[col].source) {
                    $.fn.jexcel.defaults[id].columns[col].source = [];
                }
                if (! options.columns[col].options) {
                    $.fn.jexcel.defaults[id].columns[col].options = [];
                }

                i++;
            }

            Array.prototype.push.apply(options.columns, currentColumns);
            Array.prototype.push.apply(options.colWidths, currentColWidths);
            Array.prototype.push.apply(options.colHeaders, currentColHeaders);
            Array.prototype.push.apply(options.colAlignments, currentColAlignments);
            Array.prototype.push.apply(options.colHeaderClasses, currentColHeaderClasses);

            // Visual reference
            if (direction == 1) {
                var referenceHeader = $(this).find('.jexcel-header #col-' + (columnNumber - 1));

                // Add width control
                var referenceWidthHeader = $(this).find('.jexcel-header > table > colgroup > col').eq(columnNumber);
                var referenceWidthContent = $(this).find('.jexcel-content > table > colgroup > col').eq(columnNumber);

                // Adding visual headers
                for (var col = (numOfColumns + columnNumber - 1); col >= columnNumber; col--) {
                    // Adding width control
                    $(referenceWidthHeader).after('<col width="' + options.colWidths[col] + '">')
                    $(referenceWidthContent).after('<col width="' + options.colWidths[col] + '">')
                    // Adding the header
                    var title = options.colHeaders[col];
                    var header = title ? title : $.fn.jexcel('getColumnName', col);
                    var display = (options.columns[col].type == 'hidden') ? ' style="display:none"' : '';
                    td = '<td align="' + options.colAlignments[col] + '" title="' + title + '" class="' + options.colHeaderClasses[col] + '"' + display + '>' + header + '</td>';
                    $(referenceHeader).after(td);
                }

                // Adding visual columns
                for (var row = 0; row < options.data.length; row++) {
                    // Keep the current data
                    var currentData = $.fn.jexcel.defaults[id].data[row].splice(columnNumber);

                    // Reference
                    var referenceColumn = $(this).find('tbody').find('#row-' + row).find('.c' + (columnNumber - 1));

                    for (var col = (numOfColumns + columnNumber - 1); col >= columnNumber; col--) {
                        td = $(this).jexcel('createCell', col, row);
                        $(referenceColumn).after(td);

                        // New value
                        $.fn.jexcel.defaults[id].data[row][col] = data[row] ? data[row] : '';

                        // Cell data
                        records.push({
                            col: col,
                            row: row,
                            cell: $(td),
                            newValue: $.fn.jexcel.defaults[id].data[row][col],
                            oldValue: '',
                        });
                    }

                    // Copy the data back to the main data
                    Array.prototype.push.apply($.fn.jexcel.defaults[id].data[row], currentData);
                }
            } else {
                var referenceHeader = $(this).find('.jexcel-header #col-' + columnNumber);

                // Add width control
                var referenceWidthHeader = $(this).find('.jexcel-header > table > colgroup > col').eq(columnNumber + 1);
                var referenceWidthContent = $(this).find('.jexcel-content > table > colgroup > col').eq(columnNumber + 1);

                // Adding visual headers
                for (var col = columnNumber; col < numOfColumns + columnNumber; col++) {
                    // Adding width control
                    $(referenceWidthHeader).before('<col width="' + options.colWidths[col] + '">')
                    $(referenceWidthContent).before('<col width="' + options.colWidths[col] + '">')
                    // Adding the header
                    var title = options.colHeaders[col];
                    var header = title ? title : $.fn.jexcel('getColumnName', col);
                    var display = (options.columns[col].type == 'hidden') ? ' style="display:none"' : '';
                    td = '<td align="' + options.colAlignments[col] + '" title="' + title + '" class="' + options.colHeaderClasses[col] + '"' + display + '>' + header + '</td>';
                    $(referenceHeader).before(td);
                }

                // Adding visual columns
                for (var row = 0; row < options.data.length; row++) {
                    // Keep the current data
                    var currentData = $.fn.jexcel.defaults[id].data[row].splice(columnNumber);

                    // Reference
                    var referenceColumn = $(this).find('tbody').find('#row-' + row).find('.c' + columnNumber);

                    // Add cells
                    for (var col = columnNumber; col < numOfColumns + columnNumber; col++) {
                        td = $(this).jexcel('createCell', col, row);
                        $(referenceColumn).before(td);

                        // New value
                        $.fn.jexcel.defaults[id].data[row][col] = data[row] ? data[row] : '';

                        // Cell data
                        records.push({
                            col: col,
                            row: row,
                            cell: $(td),
                            newValue: $.fn.jexcel.defaults[id].data[row][col],
                            oldValue: '',
                        });
                    }

                    // Copy the data back to the main data
                    Array.prototype.push.apply($.fn.jexcel.defaults[id].data[row], currentData);
                }
            }

            // Remove table references
            $(this).jexcel('updateTableReferences', columnNumber, 0);

            // Update selection
            $(this).jexcel('updateSelection', records[0].cell, records[records.length-1].cell);

            // Update cells
            $(this).jexcel('loadCells', records, true);

            // Keeping history of changes
            $(this).jexcel('setHistory', records, {
                type:'insertColumn',
                mixed: mixed,
                properties: properties,
                columnNumber: columnNumber,
            });

            // Events
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof(options.oninsertcolumn) == 'function') {
                    options.oninsertcolumn($(this));
                }
            }
        }
    },

    /**
     * Insert a new row
     * 
     * @param mixed - number of blank lines to be insert or a single array with the data of the new row
     * @param rowNumber
     * @return void
     */
    insertRow : function(mixed, rowNumber) {
        // Main
        var main = $(this);

        // Id
        var id = $(this).prop('id');

        // Main configuration
        var options = $.fn.jexcel.defaults[id];

        // Configuration
        if (options.allowInsertRow == true) {
            // Records
            var records = [];

            // Data to be insert
            var data = [];

            // The insert could be lead by number of rows or the array of data
            if (mixed > 0) {
                var numOfRows = mixed;
            } else {
                var numOfRows = 1;

                if (mixed) {
                    data = mixed;
                }
            }

            // Current column number
            var lastRow = options.data.length;

            // Append direction
            var direction = 0;

            // Where?
            if (rowNumber == undefined || rowNumber >= parseInt(lastRow)) {
                // Default is the end of the table
                rowNumber = parseInt(lastRow);

                // Direction
                var direction = 1;
            }

            // Keep the current data
            var currentData = $.fn.jexcel.defaults[id].data.splice(rowNumber);

            // Add row to the table body
            if (direction == 1) {
                var referenceRow = $(this).find('tbody').find('#row-' + (rowNumber - 1));
            } else {
                var referenceRow = $(this).find('tbody').find('#row-' + rowNumber);
            }

            // Creat row closure
            var createRow = function() {
                // New line of data to be append in the table
                var tr = document.createElement('tr');

                // Index column
                $(tr).append('<td class="jexcel_label"></td>');

                // New row
                $.fn.jexcel.defaults[id].data[row] = [];

                // New columns
                for (var col = 0; col < $.fn.jexcel.defaults[id].colHeaders.length; col++) {
                    // New Data
                    $.fn.jexcel.defaults[id].data[row][col] = data[col] ? data[col] : '';

                    // New column of data to be append in the line
                    td = $(main).jexcel('createCell', col, row);

                    // Add column to the row
                    $(tr).append(td);

                    // Cell data
                    records.push({
                        col: col,
                        row: row,
                        cell: $(td),
                        newValue: $.fn.jexcel.defaults[id].data[row][col],
                        oldValue: '',
                    });
                }

                return tr;
            }

            // Adding lines
            if (direction == 1) {
                for (var row = (numOfRows + rowNumber - 1); row >= rowNumber; row--) {
                    $(referenceRow).after(createRow());
                }
            } else {
                for (var row = rowNumber; row < (numOfRows + rowNumber); row++) {
                    $(referenceRow).before(createRow());
                }
            }

            // Copy the data back to the main data
            Array.prototype.push.apply($.fn.jexcel.defaults[id].data, currentData);

            // Remove table references
            $(this).jexcel('updateTableReferences');

            // Update selection
            $(this).jexcel('updateSelection', records[0].cell, records[records.length-1].cell);

            // Update cells
            $(this).jexcel('loadCells', records, true);

            // Keeping history of changes
            $(this).jexcel('setHistory', records, {
                type:'insertRow',
                mixed: mixed,
                rowNumber: rowNumber,
            });

            // Insert events
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof(options.oninsertrow) == 'function') {
                    options.oninsertrow($(this));
                }
            }
        }
    },

    /**
     * Delete a column by number
     * 
     * @param integer columnNumber - column show be excluded
     * @return void
     */
    deleteColumn : function(columnNumber, numOfColumns) {
        // Id
        var id = $(this).prop('id');

        // Main configuration
        var options = $.fn.jexcel.defaults[id];

        // Global Configuration
        if (options.allowDeleteColumn == true) {
            // Keep history
            var records = [];

            // Delete column definitions
            if (columnNumber == undefined) {
                var number = $(this).jexcel('getSelectedColumns');

                if (! number[0]) {
                    columnNumber = $.fn.jexcel.defaults[id].data[0].length - 1;
                    numOfColumns = 1;
                } else {
                    columnNumber = parseInt(number[0]);
                    numOfColumns = parseInt(number.length);
                }
            }

            // Remove the last column
            if (columnNumber == -1) {
                columnNumber = $.fn.jexcel.defaults[id].data[0].length - 1;
            }

            // Minimum of columns to be delete is 1
            if (numOfColumns != parseInt(numOfColumns)) {
                numOfColumns = 1;
            }

            // Can't delete more than the limit of the table
            if (numOfColumns > options.data[0].length - columnNumber) {
                numOfColumns = options.data[0].length - columnNumber;
            }

            // Can't remove the last column
            if (options.data[0].length > 1) {
                if (parseInt(columnNumber) > -1) {
                    // Delete the column properties
                    var columns = options.columns.splice(columnNumber, numOfColumns);
                    var colHeaders = options.colHeaders.splice(columnNumber, numOfColumns);
                    var colWidths = options.colWidths.splice(columnNumber, numOfColumns);
                    var colAlignments = options.colAlignments.splice(columnNumber, numOfColumns);
                    var colHeaderClasses = options.colHeaderClasses.splice(columnNumber, numOfColumns);

                    // Delete data from source
                    for (var row = 0; row < $.fn.jexcel.defaults[id].data.length; row++) {
                        // Keep columns data
                        data = options.data[row].splice(columnNumber, numOfColumns);

                        for (var col = 0; col < numOfColumns; col++) {
                            if (! records[col]) {
                                records[col] = [];
                            }
                            records[col].push(data[col]);
                        }
                    }

                    // Remove visual columns
                    var isColumnRemoved = false;

                    for (var col = columnNumber; col < columnNumber + numOfColumns; col++) {
                        // Remove header
                        var header =  $(this).find('thead #col-' + col);
                        if ($(header).hasClass('selected')) {
                            isColumnRemoved = true;
                        }
                        $(header).remove();

                        // Remove column
                        $(this).find('tbody .c' + col).remove();
                    }

                    // Update table references
                    $(this).jexcel('updateTableReferences', columnNumber, 0);

                    // Select first row any of the selected rows were excluded
                    if (isColumnRemoved) {
                        // Get cell objects
                        var o1 = $(this).find('#0-0');
                        var o2 = $(this).find('#0-0');

                        // Update selection
                        $(this).jexcel('updateSelection', o1, o2);
                    }
                }

                // Keeping history of changes
                $(this).jexcel('setHistory', null, {
                    type:'deleteColumn',
                    columnNumber:columnNumber,
                    numOfColumns:numOfColumns,
                    columnData:records,
                    columnProperties: {
                        columns:columns,
                        colHeaders:colHeaders,
                        colWidths:colWidths,
                        colAlignments:colAlignments,
                        colHeaderClasses:colHeaderClasses,
                    }
                });

                // Delete
                if ($.fn.jexcel.ignoreEvents != true) {
                    if (typeof(options.ondeletecolumn) == 'function') {
                        options.ondeletecolumn($(this), columnNumber, numOfColumns, records);
                    }
                }
            } else {
                console.error('It is not possible to delete the last column');
            }
        }
    },

    /**
     * Delete a row by number
     * 
     * @param integer rowNumber - row number show be excluded
     * @param integer numOfRows - number of lines
     * @return void
     */
    deleteRow : function(rowNumber, numOfRows) {
        // Id
        var id = $(this).prop('id');

        // Main configuration
        var options = $.fn.jexcel.defaults[id];

        // Global Configuration
        if (options.allowDeleteRow == true) {
            // Keep history
            var records = [];

            // Delete row definitions
            if (rowNumber == undefined) {
                var number = $(this).jexcel('getSelectedRows');

                if (! number[0]) {
                    rowNumber = $.fn.jexcel.defaults[id].data.length - 1;
                    numOfRows = 1;
                } else {
                    rowNumber = parseInt(number[0]);
                    numOfRows = parseInt(number.length);
                }
                
                rowNumber = parseInt(number[0]);
                numOfRows = parseInt(number.length);
            }

            // Remove the last column
            if (rowNumber == -1) {
                rowNumber = $.fn.jexcel.defaults[id].data.length - 1;
            }

            if (numOfRows != parseInt(numOfRows)) {
                numOfRows = 1;
            }

            // Can't delete more than the limit of the table
            if (numOfRows > options.data.length - rowNumber) {
                numOfRows = options.data.length - numOfRows;
            }

            // Can't remove the last row
            if (options.data.length > 1) {
                if (parseInt(rowNumber) > -1) {
                    // Test
                    var isRowRemoved = false;

                    // Remove rows
                    for (var row = rowNumber; row < rowNumber + numOfRows; row++) {
                        // Keep row data
                        records[row] = $.fn.jexcel.defaults[id].data[row];

                        // Remove visualy
                        var line = $(this).find('#row-' + row);

                        // The current row was excluded
                        if ($(line).hasClass('selected')) {
                            isRowRemoved = true;
                        }

                        $(line).remove();
                    }

                    // Remove from source
                    $.fn.jexcel.defaults[id].data.splice(parseInt(rowNumber), numOfRows);

                    // Remove table references
                    $(this).jexcel('updateTableReferences', 0, rowNumber);

                    // Select first row any of the selected rows were excluded
                    if (isRowRemoved) {
                        // Get cell objects
                        var o1 = $(this).find('#0-0');
                        var o2 = $(this).find('#0-0');

                        // Update selection
                        $(this).jexcel('updateSelection', o1, o2);
                    }

                    // Keeping history of changes
                    $(this).jexcel('setHistory', null, {
                        type:'deleteRow',
                        rowNumber: rowNumber,
                        numOfRows: numOfRows,
                        rowData: records
                    });

                    // Delete events
                    if ($.fn.jexcel.ignoreEvents != true) {
                        if (typeof(options.ondeleterow) == 'function') {
                            options.ondeleterow($(this), rowNumber, numOfRows, records);
                        }
                    }
                }
            } else {
                console.error('It is not possible to delete the last row');
            }
        }
    },

    /**
     * Set the column width
     * 
     * @param column - column number (first column is: 0)
     * @param width - new column width
     */
    setWidth : function (column, width) {
        var id = $(this).prop('id');

        if (width > 0) {
            // In case the column is an object
            if (typeof(column) == 'object') {
                column = $(column).prop('id').split('-');
                column = column[0];
            }

            var colWidths = $(this).find('.jexcel-header col');
            $(colWidths[parseInt(column) + 1]).prop('width', width);
            var colWidths = $(this).find('.jexcel-content col');
            var currentWidth = $(colWidths[parseInt(column) + 1]).prop('width');
            $(colWidths[parseInt(column) + 1]).prop('width', width);

            // Update config container
            $.fn.jexcel.defaults[id].colWidths[column] = width;

            // Keeping history of changes
            $(this).jexcel('setHistory', null, {
                type:'setWidth',
                column: column,
                fr: currentWidth,
                to: width,
            });
        }
    },

    /**
     * Get the column width
     * 
     * @param column - column number (first column is: 0)
     * @return width - current column width
     */
    getWidth : function (column) {
        // Id
        var id = $(this).prop('id');

        if (! column) {
            // Get all headers
            return $.fn.jexcel.defaults[id].colWidths;
        } else {
            // In case the column is an object
            if (typeof(column) == 'object') {
                column = $(column).prop('id').split('-');
                column = column[0];
            }

            return $.fn.jexcel.defaults[id].colWidths[column];
        }
    },

    /**
     * Get the column title
     * 
     * @param column -  column number (first column is: 0)
     * @param title - new column title
     */
    getHeader : function (column) {
        var col = $(this).find('.jexcel_headers #col-' + column);
        if (col.length) {
            return $(col).html();
        }
    },

    /**
     * Get the column title
     * 
     * @param column - column number (first column is: 0)
     * @param title - new column title
     */
    getHeaders : function (convertToText) {
        var col = $(this).find('.jexcel_headers > td');

        var txt = [];
        $.each(col, function(k, v) {
            if (k > 0) {
                txt.push($(v).text());
            }
        });

        if (convertToText == true) {
            txt = txt.join(',');
        }

        return txt;
    },

    /**
     * Set the column title
     * 
     * @param column - column number (first column is: 0)
     * @param title - new column title
     */
    setHeader : function (column, title) {
        var col = $(this).find('thead #col-' + column);
        if (col.length) {
            if (! title) {
                title = prompt('Column name', $(col).text())
            }
            var colHtml = $(col).html();
            colHtml = colHtml.replace($(col).text(), title);
            $(col).html(colHtml);
        }
    },

    /**
     * Update column source for dropboxes
     */
    setSource : function (column, source) {
        // In case the column is an object
        if (typeof(column) == 'object') {
            column = $(column).prop('id').split('-');
            column = column[0];
        }

        // Id
        var id = $(this).prop('id');

        // Update defaults
        $.fn.jexcel.defaults[id].columns[column].source = source;
    },

    /**
     * After change
     */
    afterChange : function(ignoreSpare) {
        // Id
        var id = $(this).prop('id');

        if (! ignoreSpare) {
            // Spare check
            $(this).jexcel('spareCheck');
        }

        // Update position
        $(this).jexcel('updateCornerPosition');

        // Update settings
        $(this).jexcel('updateSettings');

        // After Changes
        if ($.fn.jexcel.ignoreEvents != true) {
            if (typeof($.fn.jexcel.defaults[id].onafterchange) == 'function') {
                $.fn.jexcel.defaults[id].onafterchange($(this));
            }
        }
    },

    /**
     * Helper function to copy data using the corner icon
     */
    copyData : function(o, d) {
        // Get data from all selected cells
        var data = $(this).jexcel('getData', true);

        // Cells
        var px = parseInt(o[0]);
        var ux = parseInt(d[0]);
        var py = parseInt(o[1]);
        var uy = parseInt(d[1]);

        // Records
        var records = []; 

        // Copy data procedure
        var posx = 0;
        var posy = 0;
        for (j = py; j <= uy; j++) {
            // Controls
            if (data[posy] == undefined) {
                posy = 0;
            }
            posx = 0;

            // Data columns
            for (i = px; i <= ux; i++) {
                // Column
                if (data[posy] == undefined) {
                    posx = 0;
                } else if (data[posy][posx] == undefined) {
                    posx = 0;
                }

                // Get cell
                cell = $(this).find('#' + i + '-' + j);

                // Update non-readonly
                if ($(cell).length && ! $(cell).hasClass('readonly')) {
                    // Keep cells history
                    records.push({
                        col: i,
                        row: j,
                        cell: $(cell),
                        newValue: data[posy][posx],
                        oldValue: $.fn.jexcel.defaults[$.fn.jexcel.current].data[j][i],
                    });
                }
                posx++;
            }
            posy++;
        }

        // Save data
        if (records.length > 0) {
            $(this).jexcel('updateCells', records);
        }
    },

    /**
     * Sort data and reload table
     */
    orderBy :  function(column, order) {
        if (column >= 0) {
            // Identify thead container
            var c = $(this).find('#col-' + column).parent();

            // No order specified then toggle order
            if (! (order == '0' || order == '1')) {
                var d = $(c).find('.arrow-down');
                if ($(d).length > 0) {
                    order = 1;
                } else {
                    order = 0;
                }
            }

            // Remove styling
            $(c).find('.arrow-down').remove();
            $(c).find('.arrow-up').remove();
            $(c).find('td').css('text-decoration', 'none');

            // Add style on specified column
            if (order == 1) {
                $(c).find('#col-' + column).append('<span class="arrow-up"></span>').css('text-decoration', 'underline');
            } else {
                $(c).find('#col-' + column).append('<span class="arrow-down"></span>').css('text-decoration', 'underline');
            }

           // Hide corner
            $('.jexcel_corner').css('top', '-200px');
            $('.jexcel_corner').css('left', '-200px');

            // Id
            var id = $(this).prop('id');
            var options = $.fn.jexcel.defaults[id];

            Array.prototype.sortBy = function(p, o) {
                return this.slice(0).sort(function(a, b) {
                    if (Number(a[p]) == a[p] && Number(b[p]) == b[p]) {
                        var valueA = Number(a[p]);
                        var valueB = Number(b[p]);
                    } else {
                        var valueA = a[p].toLowerCase();
                        var valueB = b[p].toLowerCase();
                    }

                    if (! o) {
                        return (valueA > valueB) ? 1 : (valueA < valueB) ? -1 : 0;
                    } else {
                        return (valueA > valueB) ? -1 : (valueA < valueB) ? 1 : 0;
                    }
                });
            }

            // Reorder
            options.data = options.data.sortBy(column, order);

            // TODO: Reset data
            $(this).jexcel('setData', null, true);

            // On sort event
            if ($.fn.jexcel.ignoreEvents != true) {
                if (typeof($.fn.jexcel.defaults[id].onsort) == 'function') {
                    $.fn.jexcel.defaults[id].onsort($(this), column, order);
                }
            }

            return true;
        }
    },

    /**
     * Update formula chain
     */
    updateFormulas : function(columnId, value) {
        if (typeof excelFormulaUtilities == 'object') {
            // Id
            var id = $(this).prop('id');

            // Create chain
            var regex = /([A-Z]+[0-9]+)*/g;

            // Formula
            var formula = excelFormulaUtilities.formula2JavaScript(value);

            // Elements
            var elements = formula.match(regex);

            $.each(elements, function(k, element) {
                if (element) {
                    // Get coords reference
                    elementId = $.fn.jexcel('getIdFromColumnName', element);
                    if (! $.fn.jexcel.defaults[id].formula[elementId]) {
                        $.fn.jexcel.defaults[id].formula[elementId] = [];
                    }

                    if (elementId == columnId) {
                        console.error('Self reference error: ' + element);
                    } else {
                        $.fn.jexcel.defaults[id].formula[elementId].push(columnId);
                    }
                }
            });
        }
    },

    /**
     * Apply formula to all columns in the table
     */
    formula : function(column) {
        if (typeof excelFormulaUtilities == 'object') {
            // Id
            var id = $(this).prop('id');

            // Main
            var main = $(this);

            // Any error
            var error = false;

            // Get all formulas in the chain
            if ($.fn.jexcel.defaults[id].formula[column]) {
                // Go through all cells from the chain
                $.each($.fn.jexcel.defaults[id].formula[column], function(k, v) {
                    // Cell identification
                    var cellId = v.split('-');

                    // Get formula
                    var formula = $.fn.jexcel.defaults[id].data[cellId[1]][cellId[0]];

                    // Check if still a formula
                    if (formula.substr(0,1) != '=') {
                        // Remove any existing calculation error
                        $(main).find('#' + v).removeClass('error');
                        // Remove this column from the formula chain
                        $.fn.jexcel.defaults[id].formula[column].splice(k, 1);
                    } else {
                        $(main).jexcel('executeFormula', v);
                    }

                    // Recursive formula
                    $(main).jexcel('formula', v);
                });
            }
        } else {
            console.error('excelFormulaUtilities lib not included');
        }
    },

    /**
     * Run the formula for one given cell
     */
    executeFormula : function (cellId) {
        // Not calculate and not show error if excelFormulaUtilities is not defined
        if (typeof excelFormulaUtilities !== 'object')
            return;

        // Id
        var id = $(this).prop('id');

        // Object
        var cell = $(this).find('#' + cellId);

        // Cell identification
        var cellId = cellId.split('-');

        // Get formula
        var formula = $.fn.jexcel.defaults[id].data[cellId[1]][cellId[0]];

        // Convert formula to javascript
        try {
            var value = eval(excelFormulaUtilities.formula2JavaScript(formula));
        } catch (e) {
            var value = null;
        }

        // Set value
        if (value === null || value == undefined || (''+value) == 'NaN' || typeof value === 'function') {
            // New cell value
            value = '<input type="hidden" value="' + formula + '">#ERROR';
            // Add class error to the cell
            $(cell).addClass('error');
            // Update cell content
            $(cell).html(value);
        } else {
            // Update variables
            var letter = $.fn.jexcel('getColumnName', cellId[0])  + (parseInt(cellId[1]) + 1);
            // Update global variable
            if (value == Number(value)) {
                $.fn.jexcel.defaults[id].values[letter] = Number(value);
            } else {
                $.fn.jexcel.defaults[id].values[letter] = value;
            }
            // Update window
            window[letter] = new $.fn.jexcel.factory(letter, $.fn.jexcel.defaults[id].values[letter]);
            // New cell value
            value = '<input type="hidden" value="' + formula + '">' + value;
            // Remove any error class
            $(cell).removeClass('error');
            // Update cell content
            $(cell).html(value);
        }
    },

    /**
     * Update all cells with formulas
     * 
     * @return void
     */
    updateAllCellsWithFormulas : function() {
        if (typeof excelFormulaUtilities == 'object') {
            // Id
            var id = $(this).prop('id');

            // Columns
            var columnIds = Object.keys($.fn.jexcel.defaults[id].formula);

            // Get all formulas in the chain
            if (columnIds.length > 0) {
                for (var i = 0; i < columnIds.length; i++) {
                    $(this).jexcel('formula', columnIds[i]);
                }
            }
        }
    },

    /**
     * Multi-utility helper
     * 
     * @param object options { action: METHOD_NAME }
     * @return mixed
     */
    helper : function (options) {
        var data = [];
        if (typeof(options) == 'object') {
            // Return a empty bidimensional array
            if (options.action == 'createEmptyData') {
                var x = options.cols || 10;
                var y = options.rows || 100;
                for (j = 0; j < y; j++) {
                    data[j] = [];
                    for (i = 0; i < x; i++) {
                        data[j][i] = '';
                    }
                }
            }
        }

        return data;
    },

    /**
     * Download CSV table
     * 
     * @return null
     */
    download : function () {
        // Get table id
        var id = $(this).prop('id');
        // Increment and get the current history index
        var options = $.fn.jexcel.defaults[id];
        // Data
        var data = '';

        // Get headers if applicable
        if (options.csvHeaders == true) {
            data = $(this).jexcel('getHeaders', true) + "\n";
        }

        // Get data
        data += $(this).jexcel('copy', false, ',', true);

        // Download element
        var blob = new Blob(["\uFEFF"+data], {type: 'text/csv;charset=utf-8;'});

        // IE Compatibility
        if (window.navigator && window.navigator.msSaveOrOpenBlob) {
            window.navigator.msSaveOrOpenBlob(blob, options.csvFileName + '.csv');
        } else {
            var pom = document.createElement('a');
            var url = URL.createObjectURL(blob);
            pom.href = url;
            pom.setAttribute('download', options.csvFileName + '.csv');
            document.body.appendChild(pom);
            pom.click();
            document.body.removeChild(pom);
        }
    },

    /**
     * Initializes a new history record for undo/redo
     * 
     * @return null
     */
    setHistory : function(changes, action) {
        if ($.fn.jexcel.ignoreHistory != true) {
            var main = $(this);

            var id = $(this).prop('id');

            // Increment and get the current history index
            var index = ++$.fn.jexcel.defaults[id].historyIndex;

            // Slice the array to discard undone changes
            var history = ($.fn.jexcel.defaults[id].history = $.fn.jexcel.defaults[id].history.slice(0, index + 1));

            if (action) {
                history[index] = {
                    action:action
                }
            } else {
                // Create history slot
                history[index] = {
                    firstSelected: changes[0].cell,
                    lastSelected: changes[changes.length - 1].cell,
                    cellChanges: changes,
                };
            }
        }
    },

    /**
     * Undo last action
     */
    undo : function () {
        var id = $(this).prop('id');

        // Ignore events and history
        var ignoreEvents = $.fn.jexcel.ignoreEvents ? true : false;
        var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

        $.fn.jexcel.ignoreEvents = true;
        $.fn.jexcel.ignoreHistory = true;

        // Records
        var records = [];

        // Update cells
        if ($.fn.jexcel.defaults[id].historyIndex >= 0) {
            // History
            var historyRecord = $.fn.jexcel.defaults[id].history[$.fn.jexcel.defaults[id].historyIndex--];

            if (historyRecord.action) {
                // Remove the inserted row
                if (historyRecord.action.type == 'insertRow') {
                    var rowNumber = historyRecord.action.rowNumber;
                    var numOfRows = historyRecord.action.mixed > 0 ? historyRecord.action.mixed : 1;

                    $(this).jexcel('deleteRow', rowNumber, numOfRows);

                // Create the deleted row back
                } else if (historyRecord.action.type == 'deleteRow') {
                    var rowNumber = historyRecord.action.rowNumber;
                    var numOfRows = historyRecord.action.numOfRows;

                    for (var row = rowNumber; row < rowNumber + numOfRows; row++) {
                        $(this).jexcel('insertRow', historyRecord.action.rowData[row], row);
                    }

                // Remove the inserted column
                } else if (historyRecord.action.type == 'insertColumn') {
                    var columnNumber = historyRecord.action.columnNumber;
                    var numOfColumns = historyRecord.action.mixed > 0 ? historyRecord.action.mixed : 1;

                    $(this).jexcel('deleteColumn', columnNumber, numOfColumns);

                // Create the deleted column back
                } else if (historyRecord.action.type == 'deleteColumn') {
                    var columnNumber = historyRecord.action.columnNumber;
                    var numOfColumns = historyRecord.action.numOfColumns;
                    var columnProperties = historyRecord.action.columnProperties;

                    var columnIndex = 0;
                    for (var col = columnNumber; col < columnNumber + numOfColumns; col++) {
                        $(this).jexcel('insertColumn',
                            historyRecord.action.columnData[columnIndex], {
                                columns: [columnProperties.columns[columnIndex]],
                                colHeaders: [columnProperties.colHeaders[columnIndex]],
                                colWidths: [columnProperties.colWidths[columnIndex]],
                                colAlignments: [columnProperties.colAlignments[columnIndex]],
                                colHeaderClasses: [columnProperties.colHeaderClasses[columnIndex]],
                            }, col );
                        columnIndex++;
                    }

                // Move the row back to the original position
                } else if (historyRecord.action.type == 'moveRow') {
                    $(this).jexcel('moveRow', historyRecord.action.to, historyRecord.action.fr);

                // Change the CSS back
                } else if (historyRecord.action.type == 'setStyle') {
                    $(this).jexcel('setStyle', historyRecord.action.fr);

                // Change the width
                } else if (historyRecord.action.type == 'setWidth') {
                    $(this).jexcel('setWidth', historyRecord.action.column, historyRecord.action.fr);
                }
            } else {
                // Redo for changes in cells
                for (var i = 0; i < historyRecord.cellChanges.length; i++) {
                    // Keep cells history
                    records.push({
                        col: historyRecord.cellChanges[i].col,
                        row: historyRecord.cellChanges[i].row,
                        cell: historyRecord.cellChanges[i].cell,
                        newValue: historyRecord.cellChanges[i].oldValue,
                        oldValue: historyRecord.cellChanges[i].newValue,
                    });
                }

                // Select cell
                $.fn.jexcel.selectedCell = historyRecord.firstSelected;

                // Update selection
                $(this).jexcel('updateSelection', historyRecord.firstSelected, historyRecord.lastSelected);

                // Update cells
                $(this).jexcel('loadCells', records);
            }
        }

        $.fn.jexcel.ignoreEvents = ignoreEvents;
        $.fn.jexcel.ignoreHistory = ignoreHistory;

        $(this).jexcel('updateAllCellsWithFormulas');
        $(this).jexcel('afterChange');
    },

    /**
     * Redo previously undone action
     */
    redo : function () {
        var id = $(this).prop('id');

        var ignoreEvents = $.fn.jexcel.ignoreEvents ? true : false;
        var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

        $.fn.jexcel.ignoreEvents = true;
        $.fn.jexcel.ignoreHistory = true;

        if ($.fn.jexcel.defaults[id].historyIndex < $.fn.jexcel.defaults[id].history.length - 1) {
            // History
            var historyRecord = $.fn.jexcel.defaults[id].history[++$.fn.jexcel.defaults[id].historyIndex];

            if (historyRecord.action) {
                // Insert row
                if (historyRecord.action.type == 'insertRow') {
                    // Add the row back
                    $(this).jexcel('insertRow', historyRecord.action.mixed, historyRecord.action.rowNumber);

                // Delete row
                } else if (historyRecord.action.type == 'deleteRow') {
                    // Add the row back
                    $(this).jexcel('deleteRow', historyRecord.action.rowNumber, historyRecord.action.numOfRows);

                // Insert column
                } else if (historyRecord.action.type == 'insertColumn') {
                    // Add the row back
                    $(this).jexcel('insertColumn', historyRecord.action.mixed, historyRecord.action.properties, historyRecord.action.columnNumber);

                // Delete column
                } else if (historyRecord.action.type == 'deleteColumn') {
                    // Add the row back
                    $(this).jexcel('deleteColumn', historyRecord.action.columnNumber, historyRecord.action.numOfColumns);

                // Move row
                } else if (historyRecord.action.type == 'moveRow') {
                    // Redo the move
                    $(this).jexcel('moveRow', historyRecord.action.fr, historyRecord.action.to);

                // Redefine style
                } else if (historyRecord.action.type == 'setStyle') {
                    $(this).jexcel('setStyle', historyRecord.action.to);

                // Change the width
                } else if (historyRecord.action.type == 'setWidth') {
                    $(this).jexcel('setWidth', historyRecord.action.column, historyRecord.action.to);
                }
            } else {
                // Select cell
                $.fn.jexcel.selectedCell = historyRecord.firstSelected;

                // Update selection
                $(this).jexcel('updateSelection', historyRecord.firstSelected, historyRecord.lastSelected);

                // Update cells
                $(this).jexcel('loadCells', historyRecord.cellChanges);
            }
        }

        $.fn.jexcel.ignoreEvents = ignoreEvents;
        $.fn.jexcel.ignoreHistory = ignoreHistory;

        $(this).jexcel('updateAllCellsWithFormulas');
        $(this).jexcel('afterChange');
    },

    /**
     * Create cell
     */
    createCell : function(i, j) {
        // Get object identification
        var id = $(this).prop('id');

        // Main configuration
        var options = $.fn.jexcel.defaults[id];

        // Line properties
        align = options.colAlignments[i];
        width = options.colWidths[i];

        // Create cell and properties
        td = document.createElement('td');
        $(td).prop('width', width);
        $(td).prop('align', align);
        $(td).prop('id', i + '-' +j);
        $(td).addClass('c' + i);
        $(td).addClass('r' + j);

        // Hidden column
        if (options.columns[i].type == 'hidden') {
            $(td).css('display', 'none');
        } else if (options.columns[i].type == 'checkbox') {
            if (options.columns[i].readOnly == true) {
                $(td).html('<input type="checkbox" disabled="disabled">');
            } else {
                $(td).html('<input type="checkbox">');
            }
        } else if (options.columns[i].type == 'radio') {
            if (options.columns[i].readOnly == true) {
                $(td).html('<input type="radio" name="c' + i + '" disabled="disabled">');
            } else {
                $(td).html('<input type="radio" name="c' + i + '">');
            }
        }

        // Readonly
        if (options.columns[i].readOnly == true) {
            $(td).addClass('readonly', 'readonly');
        }

        // Wrap option
        if (options.wordWrap == true || options.columns[i].wordWrap == true) {
            $(td).css('white-space', 'pre-wrap');
        }

        // Add custom css class to column
        if (options.columns[i].cssClass) {
            $(td).addClass(options.columns[i].cssClass);
        }

        return $(td);
    },

    /**
     * Check for spare cols and rows
     */
    spareCheck : function() {
        // Id
        var id = $(this).prop('id');
        var test = false;

        // Update values
        var ignoreEvents = $.fn.jexcel.ignoreEvents ? true : false;
        var ignoreHistory = $.fn.jexcel.ignoreHistory ? true : false;

        // Disabled events and history
        $.fn.jexcel.ignoreEvents = true;
        $.fn.jexcel.ignoreHistory = true;

        // Sparerows and sparecols configuration
        if ($.fn.jexcel.defaults[id].minSpareCols > 0) {
            // Configuration to check the spare cells
            lastCol = ($.fn.jexcel.defaults[id].data[0]) ? $.fn.jexcel.defaults[id].data[0].length : 0;
            lastRow = $.fn.jexcel.defaults[id].data.length;
            checkPoint = lastCol - $.fn.jexcel.defaults[id].minSpareCols;
            if (checkPoint < 0) {
                checkPoint = 0;
            }
            // Check for non-black within the expected spare cells
            test = false;
            for (var rowNumber = 0; rowNumber < lastRow; rowNumber++) {
                for (var colNumber = checkPoint; colNumber < lastCol; colNumber++) {
                    if ($.fn.jexcel.defaults[id].data[rowNumber][colNumber]) {
                        test = true;
                    }
                }
            }
            // Spare is populated add new spare to keep it align with the configuration
            if (test) {
                $(this).jexcel('insertColumn', $.fn.jexcel.defaults[id].minSpareCols);
            }
        }

        if ($.fn.jexcel.defaults[id].minSpareRows > 0) {
            // Configuration to check the spare cells
            lastCol = ($.fn.jexcel.defaults[id].data[0]) ? $.fn.jexcel.defaults[id].data[0].length : 0;
            lastRow = $.fn.jexcel.defaults[id].data.length;
            checkPoint = lastRow - $.fn.jexcel.defaults[id].minSpareRows;
            if (checkPoint < 0) {
                checkPoint = 0;
            }
            // Check for non-black within the expected spare cells
            test = false;
            for (var rowNumber = checkPoint; rowNumber < lastRow; rowNumber++) {
                for (var colNumber = 0; colNumber < lastCol; colNumber++) {
                    if ($.fn.jexcel.defaults[id].data[rowNumber][colNumber]) {
                        test = true;
                    }
                }
            }
            // Spare is populated add new spare to keep it align with the configuration
            if (test) {
                $(this).jexcel('insertRow', $.fn.jexcel.defaults[id].minSpareCols);
            }
        }

        // Restore events and history flag
        $.fn.jexcel.ignoreEvents = ignoreEvents;
        $.fn.jexcel.ignoreHistory = ignoreHistory;
    },

    /**
     * Get seleted rows numbers
     * 
     * @return array
     */
    getSelectedRows : function() {
        var rows = [];
        // Get all selected rows
        var selectedRows = $(this).find('tbody > tr.selected');

        // Return array with all selected rows
        $.each(selectedRows, function(k, v) {
            rows.push($(v).prop('id').replace('row-', ''));
        });

        return rows;
    },

    /**
     * Get seleted rows numbers
     * 
     * @return array
     */
    getSelectedColumns : function() {
        var cols = [];
        // Get all selected rows
        var selectedColumns = $(this).find('.jexcel_headers > td.selected');

        // Return array with all selected rows
        $.each(selectedColumns, function(k, v) {
            cols.push($(v).prop('id').replace('col-', ''));
        });

        return cols;
    },

    /**
     * Get column number by class identification
     * 
     * @return array
     */
    getHeaderNumber : function(className) {
        var number = false;

        if (className) {
            // Get all selected rows
            var header = $(this).find('.jexcel_headers > td.' + className);

            if ($(header).length) {
                number = parseInt($(header).prop('id').replace('col-', ''));
            }
        }

        return number;
    },

    /**
     * Update cell references
     * 
     * @return void
     */
    updateTableReferences : function(referenceCol, referenceRow) {
        // Main
        var main = $(this);

        // Get object identification
        var id = $(this).prop('id');

        // References
        if (! referenceCol) {
            referenceCol = 0;
        }
        if (! referenceRow) {
            referenceRow = 0;
        }

        // Find cols
        if (referenceCol > -1) {
            // Get headers cells
            var headers = $(this).find('.jexcel_headers td');

            // Update all headers
            $.each(headers, function(k, v) {
                if (k > 0 && k >= referenceCol) {
                    // Update row reference
                    $(v).prop('id', 'col-' + (k - 1));

                    // Update header
                    if (! $(v).prop('title')) {
                        // Get letter
                        var header = $.fn.jexcel('getColumnName', k - 1);

                        // Update header
                        $(v).html(header)
                    }
                }
            });
        }

        // Find rows
        var rows = $(this).find('.jexcel-content tbody > tr');

        // Update all rows
        $.each(rows, function(k, v) {
            if (k >= referenceRow) {
                // Update row reference
                $(v).prop('id', 'row-' + k);

                // Find columns
                var columns = $(v).find('td');

                // Update all internal columns references
                $.each(columns, function(k1, v1) {
                    // Create new references
                    if (k1 == 0) {
                        $(v1).html(parseInt(k) + 1);
                    } else {
                        if (k1 >= referenceCol) {
                            // Remove current references
                            var coord = $(v1).prop('id').split('-');

                            // Update column reference
                            $(v1).prop('id', (k1 - 1) + '-' + k);

                            // Update row
                            $(v1).removeClass('r' + coord[1]);
                            $(v1).addClass('r' + k);

                            // Update col
                            $(v1).removeClass('c' + coord[0]);
                            $(v1).addClass('c' + (k1 - 1));

                            // Formula?
                            var val = '' + $.fn.jexcel.defaults[id].data[k][k1 - 1];

                            if (val.substr(0, 1) == '=') {
                                // Update formula references
                                if (coord[0] != k1 - 1) {
                                    // Update formulas
                                    val = $.fn.jexcel('shiftFormulaByColumn', val, k1 - 1 - coord[0]);

                                    // Update value on the cell
                                    $.fn.jexcel.defaults[id].data[k][k1 - 1] = val;
                                    $(v1).find('input').val(val);
                                }

                                if (coord[1] != k) {
                                    // Update formulas
                                    val = $.fn.jexcel('shiftFormulaByRow', val, k - coord[1]);

                                    // Update value on the cell
                                    $.fn.jexcel.defaults[id].data[k][k1 - 1] = val;
                                    $(v1).find('input').val(val);
                                }
                            }

                            // Update global values
                            var letter = $.fn.jexcel('getColumnNameFromId', [k1 - 1, k]);

                            if (val == Number(val)) {
                                $.fn.jexcel.defaults[id].values[letter] = Number(val);
                            } else {
                                $.fn.jexcel.defaults[id].values[letter] = val;
                            }

                            window[letter] = new $.fn.jexcel.factory(letter, $.fn.jexcel.defaults[id].values[letter]);
                        }
                    }
                });
            }
        });

        // Rebuild formula chain
        $.fn.jexcel.defaults[id].formula = [];

        for (var j = 0; j < $.fn.jexcel.defaults[id].data.length; j++) {
            for (var i = 0; i < $.fn.jexcel.defaults[id].data[j].length; i++) {
                var val = '' + $.fn.jexcel.defaults[id].data[j][i];
                if (val.substr(0, 1) == '=') {
                    $(main).jexcel('updateFormulas', i + '-' + j, val);
                }
            }
        }

        // Update
        $(this).jexcel('updateAllCellsWithFormulas');
        $(this).jexcel('afterChange');
    },

    /**
     * Move row
     * 
     * @return void
     */
    moveRow : function(o, d) {
        // Get object
        var id = $(this).prop('id');

        // Change data order
        $.fn.jexcel.defaults[id].data.splice(d, 0, $.fn.jexcel.defaults[id].data.splice(o, 1)[0]);

        // Clone row
        var movedRow = $(this).find('#row-' + o).clone();

        // Remove row
        $(this).find('#row-' + o).remove();

        // Move row
        if (o > d) {
            $(this).find('#row-' + d).before(movedRow);
        } else {
            $(this).find('#row-' + d).after(movedRow);
        }

        // Update references
        $(this).jexcel('updateTableReferences', 0, 0);

        // Events
        if (typeof($.fn.jexcel.defaults[id].onmoverow) == 'function') {
            $.fn.jexcel.defaults[id].onmoverow($(this), o, d);
        }

        // Keeping history of changes
        $(this).jexcel('setHistory', null, {
            type:'moveRow',
            fr: o,
            to: d,
        });
    },

    /**
     * TODO: Move column
     * 
     * @return void
     */
    moveColumn : function() {
    },

    /**
     * Get letter based on a number
     * 
     * @param integer i
     * @return string letter
     */
    getColumnName : function(i) {
        var letter = '';
        if (i > 701) {
            letter += String.fromCharCode(64 + parseInt(i / 676));
            letter += String.fromCharCode(64 + parseInt((i % 676) / 26));
        } else if (i > 25) {
            letter += String.fromCharCode(64 + parseInt(i / 26));
        }
        letter += String.fromCharCode(65 + (i % 26));

        return letter;
    },

    /**
     * Convert excel like column to jexcel id
     * 
     * @param string id
     * @return string id
     */
    getIdFromColumnName : function (id, arr) {
        // Get the letters
        var t = /^[a-zA-Z]+/.exec(id);

        if (t) {
            // Base 26 calculation
            var code = 0;
            for (var i = 0; i < t[0].length; i++) {
                code += parseInt(t[0].charCodeAt(i) - 64) * Math.pow(26, (t[0].length - 1 - i));
            }
            code--;
            // Make sure jexcel starts on zero
            if (code < 0) {
                code = 0;
            }

            // Number
            var number = parseInt(/[0-9]+$/.exec(id));
            if (number > 0) {
                number--;
            }

            if (arr == true) {
                id = [ code, number ];
            } else {
                id = code + '-' + number;
            }
        }

        return id;
    },

    /**
     * Convert jexcel id to excel like column name
     * 
     * @param string id
     * @return string id
     */
    getColumnNameFromId : function (cellId) {
        if (! Array.isArray(cellId)) {
            cellId = cellId.split('-');
        }
        return $.fn.jexcel('getColumnName', cellId[0])  + (parseInt(cellId[1]) + 1);
    },

    /**
     * Shift letters in an excel formula
     * 
     * @param string value
     * @param integer number of shifts
     */
    shiftFormulaByColumn : function(value, index) {
        // Default
        if (! index) {
            index = 1;
        }

        // Create chain
        var regex = /([A-Z]+[0-9]+)*/g;

        // Formula
        var formula = excelFormulaUtilities.formula2JavaScript(value);

        // Elements
        var elements = formula.match(regex).filter(function(n) { return n != '' }).sort(function(a, b) {
            a = $.fn.jexcel('getIdFromColumnName', a, true);
            b = $.fn.jexcel('getIdFromColumnName', b, true);

            if (index > 0) {
                return (b[0] - a[0]);
            } else {
                return (a[0] - b[0]);
            }
        });

        for (var i = 0; i < elements.length; i++) {
            if (elements[i]) {
                // Get excel-like variable
                var f = $.fn.jexcel('getIdFromColumnName', elements[i], true);
                // New letter
                var letter = f[0] + index;
                if (letter < 0) {
                    letter = 0;
                }
                // Get jexcel variable
                var t = $.fn.jexcel('getColumnName', letter);
                // Shift element from the formula
                value = value.replace(new RegExp(elements[i], "g"), t + (f[1] + 1));
            }
        }

        return value;
    },

    /**
     * Shift numbers in an excel formula
     * 
     * @param string value
     * @param integer number of shifts
     */
    shiftFormulaByRow : function(value, index) {
        // Default
        if (! index) {
            index = 1;
        }

        // Create chain
        var regex = /([A-Z]+[0-9]+)*/g;

        // Formula
        var formula = excelFormulaUtilities.formula2JavaScript(value);

        // Elements
        var elements = formula.match(regex).filter(function(n) { return n != '' }).sort(function(a, b) {
            a = $.fn.jexcel('getIdFromColumnName', a, true);
            b = $.fn.jexcel('getIdFromColumnName', b, true);

            if (index > 0) {
                return (b[1] - a[1]);
            } else {
                return (a[1] - b[1]);
            }
        });

        for (var i = 0; i < elements.length; i++) {
            if (elements[i]) {
                // Get excel-like variable
                var f = $.fn.jexcel('getIdFromColumnName', elements[i], true);
                // Get jexcel variable
                var t = '' + $.fn.jexcel('getColumnName', f[0]);
                // New number
                var number = f[1] + 1 + index;
                // New number can't be lower than 1
                if (number < 1) {
                    number = 1;
                }
                // Shift element from the formula
                value = value.replace(new RegExp(elements[i], "g"), t + number);
            }
        }

        return value;
    },


    /**
     * Get row number
     */
    row : function(cell) {
        // If blank get the current cell
        if (! cell) {
            cell = $.fn.jexcel.selectedCell;
        }

        var data = false;

        // Get object identification
        if ($(cell).length) {
            data = $(cell).prop('id').split('-');
            data = data[1];
        }

        return data;
    },

    /**
     * Get col number
     */
    col : function(cell) {
        // If blank get the current cell
        if (! cell) {
            cell = $.fn.jexcel.selectedCell
        }

        var data = false;
        // Get object identification
        if ($(cell).length) {
            data = $(cell).prop('id').split('-');
            data = data[0];
        }

        return data;
    },

    /**
     * Get events flag
     */
    getEventsFlag : function(val) {
        return $.fn.jexcel.ignoreEvents ? true : false;
    },

    /**
     * Set events flag
     */
    setEventsFlag : function(val) {
        $.fn.jexcel.ignoreEvents = val ? true : false;
    },

    /**
     * Get history flag
     */
    getHistoryFlag : function(val) {
        return $.fn.jexcel.ignoreHistory ? true : false;
    },

    /**
     * Set history flag
     */
    setHistoryFlag : function(val) {
        $.fn.jexcel.ignoreHistory = val ? true : false;
    },

    /**
     * Show filter
     */
    showFilter : function() {
        $(this).find('.jexcel-filter').css('display', '');
    },

    /**
     * Hide filter
     */
    hideFilter : function() {
        $(this).find('.jexcel-filter').css('display', 'none');
    },

    /**
     * Get cell comments
     */
    getComments : function(cell) {
        if (typeof(cell) == 'string') {
            var cell = $(this).jexcel('getCell', cell);
        }

        return $(cell).prop('title');
    },

    /**
     * Set cell comments
     */
    setComments : function(cell, comments) {
        if (typeof(cell) == 'string') {
            var cell = $(this).jexcel('getCell', cell);
        }
        $(cell).prop('title', comments);

        if (comments) {
            $(cell).addClass('jexcel_comments');
        } else {
            $(cell).removeClass('jexcel_comments');
        }
    },

    /**
     * Get config information
     */
    getConfig : function(k, v) {
        var id = $(this).prop('id');

        return $.fn.jexcel.defaults[id][k];
    },

    /**
     * Set config parameter
     */
    setConfig : function(k, v) {
        var id = $(this).prop('id');

        $.fn.jexcel.defaults[id][k] = v;
    },

    /**
     * Get number of rows'
     * 
     * @return integer
     */
    getRowCount: function() {
         // Get object
         var id = $(this).prop('id');

         // Change data order
         return $.fn.jexcel.defaults[id].data.length;
    },

    /**
     * Get meta information from cell(s)
     * 
     * @return integer
     */
    getMeta: function(cell, key) {
        // Cell
        if (! cell) {
            // Control vars
            var data = {};

            var id = $(this).prop('id');

            // Column and row length
            var x = $.fn.jexcel.defaults[id].data[0].length
            var y = $.fn.jexcel.defaults[id].data.length

            // Go through the columns to get the data
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    // Cell
                    cell = $(this).find('.jexcel-content #' + i + '-' + j);
                    // Value
                    var v = key ? $(cell).data(key) : $(cell).data();

                    // Any meta data for this column?
                    if (Object.keys(v).length) {
                        // Column name
                        var k = $.fn.jexcel('getColumnName', i)  + (parseInt(j) + 1);
                        // Value
                        data[k] = v;
                    }
                }
            }

           return data;
        } else {
            cell = $(this).jexcel('getCell', cell);

            return key ? $(cell).data(key) : $(cell).data();
        }
    },

    /**
     * Set meta information to cell(s)
     * 
     * @return integer
     */
    setMeta: function(o, k, v) {
        var main = $(this);

        if (k && v) {
            // Get object from string
            var cell = $(this).jexcel('getCell', o);
            // Set data value
            $(cell).data(k, v);
        } else {
            $.each(o, function(k, v) {
                // Get cell identifier
                var cell = Object.keys(v)[0];
                // Get column
                cell = $(main).jexcel('getCell', cell);
                // Set the data
                if (cell) {
                    $(cell).data(v[Object.keys(v)[0]]);
                }
            });
        }
    },

    /**
     * Get style information from cell(s)
     * 
     * @return integer
     */
    getStyle : function(cell, key) {
        // Cell
        if (! cell) {
            // Control vars
            var data = {};

            var id = $(this).prop('id');

            // Column and row length
            var x = $.fn.jexcel.defaults[id].data[0].length
            var y = $.fn.jexcel.defaults[id].data.length

            // Go through the columns to get the data
            for (j = 0; j < y; j++) {
                for (i = 0; i < x; i++) {
                    // Cell
                    cell = $(this).find('#' + i + '-' + j);
                    // Value
                    var v = key ? $(cell).css(key) : $(cell).attr('style');

                    // Any meta data for this column?
                    if (v) {
                        // Column name
                        var k = $.fn.jexcel('getColumnName', i)  + (parseInt(j) + 1);
                        // Value
                        data[k] = v;
                    }
                }
            }

           return data;
        } else {
            cell = $(this).jexcel('getCell', cell);

            return key ? $(cell).css(key) : $(cell).attr('style');
        }
    },

    /**
     * Set meta information to cell(s)
     * 
     * @return integer
     */
    setStyle : function(o, k, v) {
        var main = $(this);
        var styleFrom = [];
        var styleTo = [];

        if (k && v) {
            // Get object from string
            if (typeof(o) == 'string') {
                var o = $(this).jexcel('getCell', o); 
            }

            $.each(o, function(key, cell) {
                var cellId = $(cell).prop('id');

                // Current style
                var obj = [];
                var sty = $(cell).attr('style');
                obj[cellId] = '' + (sty ? sty : '');
                styleFrom.push(obj);

                // Workaround for safari ;(
                if (k == 'font-weight' && v == 'bold') {
                    if ($(cell).css(k) == '700' || $(cell).css(k) == 'bold') {
                        $(cell).css(k, '');
                    } else {
                        $(cell).css(k, v);
                    }
                } else {
                    // Change layout
                    if ($(cell).css(k) == v) {
                        $(cell).css(k, '');
                    } else {
                        $(cell).css(k, v);
                    }
                }

                // New style
                var obj = [];
                var sty = $(cell).attr('style');
                obj[cellId] = '' + (sty ? sty : '');
                styleTo.push(obj);

            });
        } else {
            $.each(o, function(k, v) {
                // Get column
                var cell = $(main).jexcel('getCell', Object.keys(v)[0]);
                // Set the data
                if (cell) {
                    // Get cell id
                    var cellId = $(cell).prop('id');

                    // Current style
                    var obj = [];
                    var sty = $(cell).attr('style');
                    obj[cellId] = '' + (sty ? sty : '');
                    styleFrom.push(obj);

                    // Change layout
                    $(cell).attr('style', v[Object.keys(v)[0]]);

                    // Current style
                    var obj = [];
                    var sty = $(cell).attr('style');
                    obj[cellId] = '' + (sty ? sty : '');
                    styleTo.push(obj);
                }
            });
        }

        // Keeping history of changes
        $(this).jexcel('setHistory', null, {
            type:'setStyle',
            fr: styleFrom,
            to: styleTo,
        });
    },

    /**
     * Default context menu
     */
    contextMenu : function(type, number, e) {
        var contextMenuContent = '';

        if (type == 'col') {
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowInsertColumn == true) {
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('insertColumn', 1, null, " + number + ")\">Insert a new column<span></span></a>";
            }
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowDeleteColumn == true) {
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('deleteColumn')\">Delete this column<span></span></a>";
            }
            contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('setHeader', " + number + ")\">Rename this column</a>";
            contextMenuContent += "<hr>";
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].columnSorting == true) {
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('orderBy', " + number + ", 0)\">Order ascending <span></span></a>";
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('orderBy', " + number + ", 1)\">Order descending <span></span></a>";
            }
            contextMenuContent += "<hr><a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('copy', true)\">Copy...<span>Ctrl + C</span></a>";
            contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('download')\">Save as...<span>Ctrl + S</span></a>";
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].about) {
                contextMenuContent += "<a onclick=\"alert('" + $.fn.jexcel.defaults[$.fn.jexcel.current].about + "')\">About<span></span></a>";
            }
        } else {
            // Default context menu for the rows
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowInsertRow == true) {
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('insertRow', 1, " + number + ")\">Insert a new row<span></span></a>";
            }
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowDeleteRow == true) {
                contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('deleteRow')\">Delete this row<span></span></a>";
            }
            contextMenuContent += "<hr>";
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].allowComments == true) {
                if (! $($.fn.jexcel.selectedCell).prop('title')) {
                    contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('setComments', $.fn.jexcel.selectedCell, prompt('Comments', ''));\">Add comments<span></span></a><hr>";
                } else {
                    contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('setComments', $.fn.jexcel.selectedCell, prompt('Comments', $($.fn.jexcel.selectedCell).prop('title')));\">Edit comments<span></span></a>";
                    contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('setComments', $.fn.jexcel.selectedCell, '');\">Reset comments<span></span></a><hr>";
                }
            }
            contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('copy', true)\">Copy...<span>Ctrl + C</span></a>";
            contextMenuContent += "<a onclick=\"jQuery('#" + $.fn.jexcel.current + "').jexcel('download')\">Save as...<span>Ctrl + S</span></a>";
            if ($.fn.jexcel.defaults[$.fn.jexcel.current].about) {
                contextMenuContent += "<a onclick=\"alert('" + $.fn.jexcel.defaults[$.fn.jexcel.current].about + "')\">About<span></span></a>";
            }
        }

        return contextMenuContent;
    },

    /**
     * Create hash based on a given string
     */
    hash : function(str) {
        var hash = 0, i, chr;

        if (str.length === 0) {
            return hash;
        } else {
            for (i = 0; i < str.length; i++) {
              chr   = str.charCodeAt(i);
              hash  = ((hash << 5) - hash) + chr;
              hash |= 0; // Convert to 32bit integer
            }
        }
        return hash;
    },

    parseTable : function(el, destroy) {
        if (el.tagName != 'TABLE') {
            console.log('Element is not a table');
        } else {
            var headers = $(el).find('thead tr');
            headers = $(headers[headers.length - 1]).find('td, th');
            options = {};
            options.columns = [];
            options.colWidths = [];
            options.colHeaders = [];
            options.colAligments = [];
            options.data = [];
            for (var i = 0; i < headers.length; i++) {
                options.columns.push({ type:'text' }),
                options.colWidths.push($(headers[i]).width() + 50),
                options.colHeaders.push($(headers[i]).text()),
                options.colAligments.push($(headers[i]).css('text-align'));
            }

            var content = $(el).find('table > tr, tbody tr');
            for (var j = 0; j < content.length; j++) {
                var row = $(content[j]).find('td');
                options.data[j] = [];
                for (var i = 0; i < row.length; i++) {
                    var value = $(row[i]).text();
                    options.data[j].push(value);
                }
            }

            var pattern = [];
            for (var i = 0; i < options.columns.length; i++) {
                test = true;
                testCalendar = true;
                pattern[i] = [];
                for (var j = 0; j < options.data.length; j++) {
                    var value = options.data[j][i];
                    if (! pattern[i][value]) {
                        pattern[i][value] = 0;
                    }
                    pattern[i][value]++;
                    if (value.length > 25) {
                        test = false;
                    }
                    if (value.length == 10) {
                        if (! (value.substr(4,1) == '-' && value.substr(7,1) == '-')) {
                            testCalendar = false;
                        }
                    } else {
                        testCalendar = false;
                    }
                }

                if (testCalendar) {
                    options.columns[i].type = 'calendar';
                } else if (test == true && Object.keys(pattern[i]).length <= parseInt(options.data.length * 0.1)) {
                    options.columns[i].type = 'dropdown';
                    options.columns[i].source = Object.keys(pattern[i]);
                }
            }

            if (destroy) {
                $(el).remove();
            }

            return options;
        }
    },

    /**
     * Destroy
     * 
     * @return array
     */
    destroy : function () {
        var id = $(this).prop('id');

        if ($.fn.jexcel.defaults[id]) {
            // Delete source
            delete $.fn.jexcel.defaults[id];

            // Remove HTML elements
            $(this).html('');

            // Remove HTML controls
            $('.jexcel_corner').remove();
            $('.jexcel_textarea').remove();
            $('.jexcel_contextmenu').remove();
            $('.jexcel_about').remove();

            // Remove selection
            if ($.fn.jexcel.current == id) {
                $.fn.jexcel.current = null;
            }
        }

        // If no other spreadsheet in the screen :: remove all elements
        if (! Object.keys($.fn.jexcel.defaults).length) {
            // Unbind events
            $(document).off('dragstart', $.fn.jexcel.dragStartControls);
            $(document).off("contextmenu", $.fn.jexcel.contextMenuControls);
            $(document).off('mousewheel', $.fn.jexcel.mouseWheelControls);
            $(document).off('mousedown touchstart', $.fn.jexcel.mouseDownControls);
            $(document).off('mouseup', $.fn.jexcel.mouseUpControls);
            $(document).off('dblclick touchend', $.fn.jexcel.doubleClickControls);
            $(document).off('mousemove', $.fn.jexcel.mouseMoveControls);
            $(document).off('mouseover', $.fn.jexcel.mouseOverControls);
            $(document).off('paste', $.fn.jexcel.pasteControls);
            $(document).off('keydown', $.fn.jexcel.keyDownControls);

            // Remove other objects
            $.fn.jexcel.defaults = null;
            $.fn.jexcel.current = null;
            $.fn.jexcel.factory = null;
            $.fn.jexcel.hash = null;
            $.fn.jexcel.style = null;
            $.fn.jexcel.selectedCorner = null;
            $.fn.jexcel.selectedHeader = null;
            $.fn.jexcel.resizeColumn = null;
            $.fn.jexcel.dragStartControls = null;
            $.fn.jexcel.contextMenuControls = null;
            $.fn.jexcel.mouseWheelControls = null;
            $.fn.jexcel.mouseDownControls = null;
            $.fn.jexcel.mouseUpControls = null;
            $.fn.jexcel.doubleClickControls = null;
            $.fn.jexcel.mouseMoveControls = null;
            $.fn.jexcel.mouseOverControls = null;
            $.fn.jexcel.pasteControls = null;
            $.fn.jexcel.keyDownControls = null;
            $.fn.jexcel.touchControls = null;
        }
    } 
};

$.fn.jexcel = function( method ) {
    if ( methods[method] ) {
        return methods[ method ].apply( this, Array.prototype.slice.call( arguments, 1 ));
    } else if ( typeof method === 'object' || ! method ) {
        return methods.init.apply( this, arguments );
    } else {
        $.error( 'Method ' +  method + ' does not exist on jQuery.tooltip' );
    }
};

})( jQuery );
