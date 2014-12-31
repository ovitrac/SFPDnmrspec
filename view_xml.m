function view_xml(xmlFileName)


pin_icon = [matlabroot, '/toolbox/matlab/icons/pin_icon.gif'];
att_icon = [matlabroot, '/toolbox/matlab/icons/greencircleicon.gif'];

fig = figure(...
    'Name','xml TableViewer', ...
    'Units', 'normalized', ...
    'NumberTitle','off', ...
    'Resize','on', ...
    'MenuBar','none', ...
    'Tag', 'xmlFigure',...
    'Visible','on');
scroll_wip_panel= uipanel('Units','normalized','Position',[0.4 0 0.6 0.7]);

table_h = uitable(scroll_wip_panel,'Units','normalized','Position',[0 0 1 0.9],'CellSelectionCallback',@cell_select,...
    'ColumnName',{'VAR_NAME','Path'});
uicontrol('Parent', scroll_wip_panel,'Units','normalized','Position',[0.2 0.92 0.6 0.06],'String','delete rows','callback',{@del_rows,table_h});



%table_h = createTable(scroll_wip_panel,{'Var name','path'},DAT,0);
search_panel= uipanel('Units','normalized','Position',[0.4 0.7 0.6 0.3],'BackgroundColor',[1,1,1]);
set(fig,'CloseRequestFcn',{@close_fun,table_h});
uicontrol('parent',search_panel,'Style','text','Units','normalized','Position',[0 0.6 0.15 0.2],'String','search:','Fontsize',12,'BackgroundColor',[1,1,1]);
%xmlFileName=[pwd,'\guitars.xml'];
%xmlFileName=[pwd,'\PVT0511011__17_11_2009_22_55_04.xml'];
h=actxserver( 'Microsoft.XMLDOM');
st=h.load(xmlFileName);
childs=h.childNodes;
kk=0;
%childs.item(0).nodeTypeString;
while ~strcmp(childs.item(kk).nodeTypeString,'element')
    kk=kk+1;
end
root=childs.item(kk);

%root.nodeTypeString
tree_root = uitreenode('v0',root.baseName, root.baseName, [], false);
drawnow;
tree_h = uitree('v0',fig,'Root', tree_root);
drawnow;
tree_h.getTree.setExpandsSelectedPaths(true);
uicontrol('parent',search_panel,'Style','edit','Units','normalized','Position',[0.155 0.6 0.4 0.2],'Callback',{@find_in_tree3,tree_h});
uicontrol('parent',search_panel,'String','ADD','Units','normalized','Position',[0.77 0.6 0.1 0.2],'Callback',{@add_to_table2,tree_h,table_h});
uicontrol('parent',search_panel,'String','NEXT','Units','normalized','Enable','off','Position',[0.56 0.6 0.2 0.2],'Callback',{@search_next,tree_h},'Tag','next');
uicontrol('parent',search_panel,'Style','text','Units','normalized','Position',[0.0 0.2 0.42 0.2],'String','parent node identifier:','Fontsize',11,'BackgroundColor',[1,1,1]);
uicontrol('parent',search_panel,'Style','edit','Units','normalized','Position',[0.45 0.2 0.5 0.2],'String','name|title','Fontsize',12,'Tag','str_handle');
set(tree_h,'DndEnabled',false);
set(tree_h,'MultipleSelectionEnabled',true);
%set(tree_h,'position',[0,300,150,100]);
set(tree_h, 'Units', 'normalized', 'position', [0 0 0.4 1.0]);
set(fig,'ResizeFcn',{@fig_resize,table_h});
%tree_h.ScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED;
child=root.childNodes;
get_child_nodes(child,tree_root);
tree_h.expand(tree_root);
%set(fig,'Visible','on');

set(fig,'ResizeFcn',{@fig_resize,[table_h]});
set(table_h,'ColumnEditable',true(1,2));
set(fig,'Position', get(fig,'Position')+0.001);
%make_sortable(table_h);
uiwait
    function get_child_nodes(child_list,parent_treenode)
        try
            % global tree_h pin_icon
            l=child_list.length;
            for ii=1:l
                this_node=child_list.item(ii-1);
                if this_node.hasChildNodes
                    nodes(ii)= uitreenode('v0',this_node.baseName, this_node.baseName, [], ~(this_node.hasChildNodes || ~isempty(this_node.attributes)));
                    drawnow;
                    get_child_nodes(this_node.childNodes,nodes(ii))
                else
                    switch this_node.nodeTypeString
                        case 'element'
                            nodes(ii)= uitreenode(this_node.nodeValue, this_node.nodeName, pin_icon,  ~(this_node.hasChildNodes || ~isempty(this_node.attributes)));
                            
                        otherwise
                            nodes(ii)= uitreenode('v0',this_node.nodeValue, this_node.nodeValue, pin_icon, ~(this_node.hasChildNodes || ~isempty(this_node.attributes)));
                    end
                end
                % nodes(ii).UserObject=this_node.nodeTypeString;
                set(nodes(ii),'UserObject',this_node.nodeTypeString);
                tree_h.add(parent_treenode, nodes);
                att=this_node.attributes;
                if ~isempty(att) && att.length>0
                    for jj=1:att.length
                        
                        att_string=att.item(jj-1).xml;
                        att_nodes(jj)=uitreenode('v0',att_string, att_string, att_icon, true);
                        %att_nodes(jj).UserObject=att.item(jj-1).nodeTypeString;
                        set(att_nodes(jj),'UserObject',att.item(jj-1).nodeTypeString);
                        set(att_nodes(jj),'UserData',att.item(jj-1).nodeName);
                        %set(att_nodes(jj),'UserData','att');
                    end
                    % try
                    tree_h.add(nodes(ii), att_nodes);
                    %nodes(ii)
                    % catch
                    %  att
                    %end
                    
                end
            end
        catch
            a=1;
        end
    end
end



function close_fun(s,e,table_h)
try
    dat=cell(table_h.Data);
    dat(all(cellfun('isempty',dat),2),:)=[];
    if isempty(dat{1,1})
        delete(gcbf);
        return
    else
        bad_names=check_names(dat(:,1));
        if ~isempty(bad_names)
            warndlg([bad_names(:)] ,'bad names');
            return
        end
    end
    path_strings=merge_str_cells(dat(:,1), repmat({'='},size(dat,1),1), dat(:,2));
    main_fig=findobj('type','figure','name','CDM EXPLORER');
    setappdata(main_fig,'path_strings',path_strings);
    delete(gcbf);
    %close(gcbf);
catch
    delete(gcbf);
end

end


function find_in_tree3(s,e,tree)

%found_nodes(1)={};
str=get(s,'String');
root_node=tree.getRoot;
childs=root_node.children;
count=0;
find_str_in_child(childs)
if exist('found_nodes','var')
    found_nodes=flipud(found_nodes);
    tree.setSelectedNode(found_nodes(end));
    drawnow;
    sel_row=tree.getTree.getLeadSelectionRow;
    sel=max(1,sel_row+0);
    tree.getTree.scrollRowToVisible(sel);
    tree.Tree.grabFocus;
    
    setappdata(gcf,'found_nodes',found_nodes);
    setappdata(gcf,'current',length(found_nodes));
    obj=findobj(gcf,'Tag','next');
    set(obj,'Enable','on');
end
%path_=get_nodes_path(found_nodes(1))
    function find_str_in_child(in_childs)
        while in_childs.hasMoreElements
            node= in_childs.nextElement;
            nodename=char(node.getName);
            if regexpi(nodename,str)
                count=count+1;
                % tree.expand(node.getParent)
                %tree.setSelectedNode(node);
                found_nodes(count)=node;
            elseif node.getAllowsChildren
                find_str_in_child(node.children);
            end
        end
    end
end

function add_to_table2(s,e,tree_h,table_h)
sel_nodes=tree_h.SelectedNodes;
l=length(sel_nodes);
new_dat=cell(l,2);
for ii=1:l
    [var_name,var_path]=get_nodes_path(sel_nodes(ii),ii);
    new_dat(ii,:)={var_name,var_path};
end
Dat=get(table_h,'Data');
Dat=[Dat;new_dat];
set(table_h,'Data',Dat);
end

function [var_name,var_path]=get_nodes_path(node,ind)
typ= get(node,'UserObject');
npath=node.getPath;


try
    switch typ
        case 'text'
            first_part=chain_nodes(npath(1:end-2));
            if length(npath)>2
                obj=findobj(gcf,'Tag','str_handle');
                str=get(obj,'String');
                ident_node= search_par(npath(end-1),str);
                if isempty(ident_node)
                    ident_node=npath(end-2).getFirstChild;
                end
                ident_name=char(ident_node.getName);
                ident_val=strtrim(char(ident_node.getFirstLeaf.getName));
                last_part=['[',ident_name,'="',ident_val,'"]/',char(npath(end-1).getName)];
                var_name=ident_val;
            else
                last_part=['/',char(npath(end).getName)];
            end
        case 'attribute'
            first_part=chain_nodes(npath(1:end-1));
            if length(npath)>1
                obj=findobj(gcf,'Tag','str_handle');
                str=get(obj,'String');
                ident_node= search_par(npath(end),str);
                if isempty(ident_node)
                    ident_node=npath(end-2).getFirstChild;
                end
                ident_name=char(ident_node.getName);
                ident_val=char(ident_node.getFirstLeaf.getName);
                
                %var_name=get(node,'UserData');
                full_var=char(node.getName);
                var_name=regexpi(full_var,'.+(?==)','match');
                var_name=var_name{:};
                last_part= ['[',ident_name,'="',ident_val,'"]','/@',var_name];
            else
                last_part=['/@',var_name];
            end
        otherwise
            first_part=chain_nodes(npath(1:end-1));
            last_part=['/',char(npath(end).getName)];
            var_name=['enter_name_',num2str(ind)];
    end
    var_path=[first_part,last_part];
catch
    switch typ
        case 'text'
            first_part=chain_nodes(npath(1:end-2));
            last_part=['/',char(npath(end-1).getName)];
            var_name=['enter_name_',num2str(ind)];
        case 'attribute'
            first_part=chain_nodes(npath(1:end-2));
            var_name=get(node,'UserData');
            last_part=['/@',nam];
        otherwise
            first_part=chain_nodes(npath(1:end-2));
            last_part=['/',char(npath(end).getName)];
            var_name=['enter_name_',num2str(ind)];
    end
    %first_part=chain_nodes(npath(1:end-1));
    var_path=[first_part,last_part];
end

end


function res=chain_nodes(npath)
res=[];
for ii=1:length(npath)-1
    res=[res,char(npath(ii).getName),'/'];
end
res=[res,char(npath(end).getName)];
end




function res= search_par(par,find_str)
%par=npath(end-1);
%find_str='name|title';


NS=par;
while 1
    NS=NS.getNextSibling;
    if isempty(NS)
        break
    else
        if regexpi(char(NS.getName),find_str)
            res=NS;
            return
        end
    end
end
NS=par;
while 1
    NS=NS.getPreviousSibling;
    if isempty(NS)
        break
    else
        if regexpi(char(NS.getName),find_str)
            res=NS;
            return
        end
    end
end
if regexpi(char(par.getName),find_str)
    res=par;
    return
end
res=[];
end

function search_next(s,e,tree)
found_nodes=getappdata(gcf,'found_nodes');
current=  getappdata(gcf,'current');
current=current-1;
if current==0
    current=length(found_nodes);
end
setappdata(gcf,'current',current);
tree.setSelectedNode(found_nodes(current));
drawnow
sel_row=tree.getTree.getLeadSelectionRow;
sel=max(1,sel_row+0);
tree.getTree.scrollRowToVisible(sel);
tree.Tree.grabFocus;

end
%('AtpDB/AtpDB[name="BiasVoltage"]/value')

function fig_resize(s,e,Tables)
pause(0.1)
for Table=Tables
    set(Table,'Units','pixels');
    t=get(Table,'Position');
    W=t(3)-32;
    n=size(get(Table,'Data'),2);
    if n==0
        n=length(get(Table,'ColumnName'));
    end
    c=cell(1,n);
    [c{:}]=deal(W/n);
    set(Table,'ColumnWidth',c);
    %pause(0.2)
    set(Table,'Units','normalized');
end
end

function cell_select(s,e)
%tabledata.Indices=e.Indices;
tabledata.selected_rows=(e.Indices(:,1));
tabledata.selected_cols=(e.Indices(:,2));
set(s,'UserData',tabledata);
end