<!DOCTYPE html>
<html class="" lang="en">
<head prefix="og: http://ogp.me/ns#">
<meta charset="utf-8">
<meta content="IE=edge" http-equiv="X-UA-Compatible">
<meta content="object" property="og:type">
<meta content="GitLab" property="og:site_name">
<meta content="README.md · master · Andre Hollstein / S2MSI" property="og:title">
<meta content="Read and use Sentinel-2 images." property="og:description">
<meta content="https://git.gfz-potsdam.de/assets/gitlab_logo-7ae504fe4f68fdebb3c2034e36621930cd36ea87924c11ff65dbcb8ed50dca58.png" property="og:image">
<meta content="https://git.gfz-potsdam.de/hollstei/S2MSI/blob/master/README.md" property="og:url">
<meta content="summary" property="twitter:card">
<meta content="README.md · master · Andre Hollstein / S2MSI" property="twitter:title">
<meta content="Read and use Sentinel-2 images." property="twitter:description">
<meta content="https://git.gfz-potsdam.de/assets/gitlab_logo-7ae504fe4f68fdebb3c2034e36621930cd36ea87924c11ff65dbcb8ed50dca58.png" property="twitter:image">

<title>README.md · master · Andre Hollstein / S2MSI · GitLab</title>
<meta content="Read and use Sentinel-2 images." name="description">
<link rel="shortcut icon" type="image/x-icon" href="/assets/favicon-075eba76312e8421991a0c1f89a89ee81678bcde72319dd3e8047e2a47cd3a42.ico" />
<link rel="stylesheet" media="all" href="/assets/application-cc91b74724180405f2eb58b823b7a68836ace4c9338b3afcd9da692fe3d9175d.css" />
<link rel="stylesheet" media="print" href="/assets/print-68eed6d8135d858318821e790e25da27b2b4b9b8dbb1993fa6765d8e2e3e16ee.css" />
<script src="/assets/application-fc596e7bbc989b689b94960672c0c4af5cb7a609cd2a9ac4e9b34bcf27f20456.js"></script>
<meta name="csrf-param" content="authenticity_token" />
<meta name="csrf-token" content="pWA3ylBV9vyv1+3lI/x+CsEEdHXeH9RVMhW81BpxJ3FMtBBrNprwK1hUxFlwlR2xgvTfdIkelY9Xl/Bs+TkiFQ==" />
<meta content="origin-when-cross-origin" name="referrer">
<meta content="width=device-width, initial-scale=1, maximum-scale=1" name="viewport">
<meta content="#474D57" name="theme-color">
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-5a9cee0e8a51212e70b90c87c12f382c428870c0ff67d1eb034d884b78d2dae7.png" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-a6eec6aeb9da138e507593b464fdac213047e49d3093fc30e90d9a995df83ba3.png" sizes="76x76" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-retina-72e2aadf86513a56e050e7f0f2355deaa19cc17ed97bbe5147847f2748e5a3e3.png" sizes="120x120" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-retina-8ebe416f5313483d9c1bc772b5bbe03ecad52a54eba443e5215a22caed2a16a2.png" sizes="152x152" />
<link color="rgb(226, 67, 41)" href="/assets/logo-d36b5212042cebc89b96df4bf6ac24e43db316143e89926c0db839ff694d2de4.svg" rel="mask-icon">
<meta content="/assets/msapplication-tile-1196ec67452f618d39cdd85e2e3a542f76574c071051ae7effbfde01710eb17d.png" name="msapplication-TileImage">
<meta content="#30353E" name="msapplication-TileColor">




<style>
  [data-user-is] {
    display: none !important;
  }
  
  [data-user-is="76"] {
    display: block !important;
  }
  
  [data-user-is="76"][data-display="inline"] {
    display: inline !important;
  }
  
  [data-user-is-not] {
    display: block !important;
  }
  
  [data-user-is-not][data-display="inline"] {
    display: inline !important;
  }
  
  [data-user-is-not="76"] {
    display: none !important;
  }
</style>

</head>

<body class="ui_charcoal" data-group="" data-page="projects:blob:show" data-project="S2MSI">
<script>
//<![CDATA[
window.gon={};gon.api_version="v3";gon.default_avatar_url="https:\/\/git.gfz-potsdam.de\/assets\/no_avatar-849f9c04a3a0d0cea2424ae97b27447dc64a7dbfae83c036c45b403392f0e8ba.png";gon.max_file_size=10;gon.relative_url_root="";gon.shortcuts_path="\/help\/shortcuts";gon.user_color_scheme="monokai";gon.award_menu_url="\/emojis";gon.current_user_id=76;
//]]>
</script>
<script>
  window.project_uploads_path = "/hollstei/S2MSI/uploads";
  window.preview_markdown_path = "/hollstei/S2MSI/preview_markdown";
</script>

<header class="navbar navbar-fixed-top navbar-gitlab with-horizontal-nav">
<div class="container-fluid">
<div class="header-content">
<button aria-label="Toggle global navigation" class="side-nav-toggle" type="button">
<span class="sr-only">Toggle navigation</span>
<i class="fa fa-bars"></i>
</button>
<button class="navbar-toggle" type="button">
<span class="sr-only">Toggle navigation</span>
<i class="fa fa-ellipsis-v"></i>
</button>
<div class="navbar-collapse collapse">
<ul class="nav navbar-nav">
<li class="hidden-sm hidden-xs">
<div class="has-location-badge search search-form">
<form class="navbar-form" action="/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><div class="search-input-container">
<div class="location-badge">This project</div>
<div class="search-input-wrap">
<div class="dropdown" data-url="/search/autocomplete">
<input type="search" name="search" id="search" placeholder="Search" class="search-input dropdown-menu-toggle" spellcheck="false" tabindex="1" autocomplete="off" data-toggle="dropdown" />
<div class="dropdown-menu dropdown-select">
<div class="dropdown-content"><ul>
<li>
<a class="is-focused dropdown-menu-empty-link">
Loading...
</a>
</li>
</ul>
</div><div class="dropdown-loading"><i class="fa fa-spinner fa-spin"></i></div>
</div>
<i class="search-icon"></i>
<i class="clear-icon js-clear-input"></i>
</div>
</div>
</div>
<input type="hidden" name="group_id" id="group_id" />
<input type="hidden" name="project_id" id="search_project_id" value="321" />
<input type="hidden" name="search_code" id="search_code" value="true" />
<script>
  gl.projectOptions = gl.projectOptions || {};
  gl.projectOptions["S2MSI"] = {
    issuesPath: "/hollstei/S2MSI/issues",
    mrPath: "/hollstei/S2MSI/merge_requests",
    name: "S2MSI"
  };
</script>
<script>
  gl.dashboardOptions = {
    issuesPath: "https://git.gfz-potsdam.de/dashboard/issues",
    mrPath: "https://git.gfz-potsdam.de/dashboard/merge_requests"
  };
</script>
<input type="hidden" name="repository_ref" id="repository_ref" value="master" />

<div class="search-autocomplete-opts hide" data-autocomplete-path="/search/autocomplete" data-autocomplete-project-id="321" data-autocomplete-project-ref="master"></div>
</form></div>

</li>
<li class="visible-sm visible-xs">
<a title="Search" aria-label="Search" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/search"><i class="fa fa-search"></i>
</a></li>
<li>
<a title="Todos" aria-label="Todos" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/dashboard/todos"><i class="fa fa-bell fa-fw"></i>
<span class="badge hidden todos-pending-count">
0
</span>
</a></li>
<li>
<a title="New project" aria-label="New project" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/projects/new"><i class="fa fa-plus fa-fw"></i>
</a></li>
<li class="header-user dropdown">
<a class="header-user-dropdown-toggle" data-toggle="dropdown" href="/u/danschef"><img width="26" height="26" class="header-user-avatar" src="https://secure.gravatar.com/avatar/7e89f8a658f61ea3b510ff308d65e0bb?s=52&amp;d=identicon" alt="7e89f8a658f61ea3b510ff308d65e0bb?s=52&amp;d=identicon" />
<span class="caret"></span>
</a><div class="dropdown-menu-nav dropdown-menu-align-right">
<ul>
<li>
<a class="profile-link" aria-label="Profile" data-user="danschef" href="/u/danschef">Profile</a>
</li>
<li>
<a aria-label="Profile Settings" href="/profile">Profile Settings</a>
</li>
<li class="divider"></li>
<li>
<a class="sign-out-link" aria-label="Sign out" rel="nofollow" data-method="delete" href="/users/sign_out">Sign out</a>
</li>
</ul>
</div>
</li>
</ul>
</div>
<h1 class="title"><a href="/u/hollstei">Andre Hollstein</a> / <a class="project-item-select-holder" href="/hollstei/S2MSI">S2MSI</a><button name="button" type="button" class="dropdown-toggle-caret js-projects-dropdown-toggle" aria-label="Toggle switch project dropdown" data-target=".js-dropdown-menu-projects" data-toggle="dropdown"><i class="fa fa-chevron-down"></i></button></h1>
<div class="header-logo">
<a class="home" title="Dashboard" id="logo" href="/"><svg width="36" height="36" class="tanuki-logo">
  <path class="tanuki-shape tanuki-left-ear" fill="#e24329" d="M2 14l9.38 9v-9l-4-12.28c-.205-.632-1.176-.632-1.38 0z"/>
  <path class="tanuki-shape tanuki-right-ear" fill="#e24329" d="M34 14l-9.38 9v-9l4-12.28c.205-.632 1.176-.632 1.38 0z"/>
  <path class="tanuki-shape tanuki-nose" fill="#e24329" d="M18,34.38 3,14 33,14 Z"/>
  <path class="tanuki-shape tanuki-left-eye" fill="#fc6d26" d="M18,34.38 11.38,14 2,14 6,25Z"/>
  <path class="tanuki-shape tanuki-right-eye" fill="#fc6d26" d="M18,34.38 24.62,14 34,14 30,25Z"/>
  <path class="tanuki-shape tanuki-left-cheek" fill="#fca326" d="M2 14L.1 20.16c-.18.565 0 1.2.5 1.56l17.42 12.66z"/>
  <path class="tanuki-shape tanuki-right-cheek" fill="#fca326" d="M34 14l1.9 6.16c.18.565 0 1.2-.5 1.56L18 34.38z"/>
</svg>

</a></div>
<div class="js-dropdown-menu-projects">
<div class="dropdown-menu dropdown-select dropdown-menu-projects">
<div class="dropdown-title"><span>Go to a project</span><button class="dropdown-title-button dropdown-menu-close" aria-label="Close" type="button"><i class="fa fa-times dropdown-menu-close-icon"></i></button></div>
<div class="dropdown-input"><input type="search" id="" class="dropdown-input-field" placeholder="Search your projects" autocomplete="off" /><i class="fa fa-search dropdown-input-search"></i><i role="button" class="fa fa-times dropdown-input-clear js-dropdown-input-clear"></i></div>
<div class="dropdown-content"></div>
<div class="dropdown-loading"><i class="fa fa-spinner fa-spin"></i></div>
</div>
</div>

</div>
</div>
</header>

<script>
  var findFileURL = "/hollstei/S2MSI/find_file/master";
</script>

<div class="page-with-sidebar">
<div class="sidebar-wrapper nicescroll">
<div class="sidebar-action-buttons">
<a class="nav-header-btn toggle-nav-collapse" title="Open/Close" href="#"><span class="sr-only">Toggle navigation</span>
<i class="fa fa-bars"></i>
</a><a class="nav-header-btn pin-nav-btn has-tooltip  js-nav-pin" title="Pin Navigation" data-placement="right" data-container="body" href="#"><span class="sr-only">Toggle navigation pinning</span>
<i class="fa fa-fw fa-thumb-tack"></i>
</a></div>
<ul class="nav nav-sidebar">
<li class="active home"><a title="Projects" class="dashboard-shortcuts-projects" href="/dashboard/projects"><span>
Projects
</span>
</a></li><li class=""><a title="Todos" href="/dashboard/todos"><span>
Todos
<span class="count">0</span>
</span>
</a></li><li class=""><a class="dashboard-shortcuts-activity" title="Activity" href="/dashboard/activity"><span>
Activity
</span>
</a></li><li class=""><a title="Groups" href="/dashboard/groups"><span>
Groups
</span>
</a></li><li class=""><a title="Milestones" href="/dashboard/milestones"><span>
Milestones
</span>
</a></li><li class=""><a title="Issues" class="dashboard-shortcuts-issues" href="/dashboard/issues?assignee_id=76"><span>
Issues
<span class="count">0</span>
</span>
</a></li><li class=""><a title="Merge Requests" class="dashboard-shortcuts-merge_requests" href="/dashboard/merge_requests?assignee_id=76"><span>
Merge Requests
<span class="count">0</span>
</span>
</a></li><li class=""><a title="Snippets" href="/dashboard/snippets"><span>
Snippets
</span>
</a></li><li class=""><a title="Help" href="/help"><span>
Help
</span>
</a></li><li class=""><a title="Profile Settings" data-placement="bottom" href="/profile"><span>
Profile Settings
</span>
</a></li></ul>

</div>
<div class="layout-nav">
<div class="container-fluid">
<div class="controls">
<div class="dropdown project-settings-dropdown">
<a class="dropdown-new btn btn-default" data-toggle="dropdown" href="#" id="project-settings-button">
<i class="fa fa-cog"></i>
<i class="fa fa-caret-down"></i>
</a>
<ul class="dropdown-menu dropdown-menu-align-right">
<li class=""><a title="Members" class="team-tab tab" href="/hollstei/S2MSI/project_members"><span>
Members
</span>
</a></li>
<li class="divider"></li>
<li>
<a data-confirm="Are you sure you want to leave the &quot;Andre Hollstein / S2MSI&quot; project?" title="Leave project" rel="nofollow" data-method="delete" href="/hollstei/S2MSI/project_members/leave">Leave Project
</a></li>
</ul>
</div>
</div>
<div class="nav-control scrolling-tabs-container">
<div class="fade-left">
<i class="fa fa-angle-left"></i>
</div>
<div class="fade-right">
<i class="fa fa-angle-right"></i>
</div>
<ul class="nav-links scrolling-tabs">
<li class="home"><a title="Project" class="shortcuts-project" href="/hollstei/S2MSI"><span>
Project
</span>
</a></li><li class=""><a title="Activity" class="shortcuts-project-activity" href="/hollstei/S2MSI/activity"><span>
Activity
</span>
</a></li><li class="active"><a title="Repository" class="shortcuts-tree" href="/hollstei/S2MSI/tree/master"><span>
Repository
</span>
</a></li><li class=""><a title="Pipelines" class="shortcuts-pipelines" href="/hollstei/S2MSI/pipelines"><span>
Pipelines
</span>
</a></li><li class=""><a title="Graphs" class="shortcuts-graphs" href="/hollstei/S2MSI/graphs/master"><span>
Graphs
</span>
</a></li><li class=""><a title="Issues" class="shortcuts-issues" href="/hollstei/S2MSI/issues"><span>
Issues
<span class="badge count issue_counter">0</span>
</span>
</a></li><li class=""><a title="Merge Requests" class="shortcuts-merge_requests" href="/hollstei/S2MSI/merge_requests"><span>
Merge Requests
<span class="badge count merge_counter">0</span>
</span>
</a></li><li class=""><a title="Wiki" class="shortcuts-wiki" href="/hollstei/S2MSI/wikis/home"><span>
Wiki
</span>
</a></li><li class="hidden">
<a title="Network" class="shortcuts-network" href="/hollstei/S2MSI/network/master">Network
</a></li>
<li class="hidden">
<a class="shortcuts-new-issue" href="/hollstei/S2MSI/issues/new">Create a new issue
</a></li>
<li class="hidden">
<a title="Builds" class="shortcuts-builds" href="/hollstei/S2MSI/builds">Builds
</a></li>
<li class="hidden">
<a title="Commits" class="shortcuts-commits" href="/hollstei/S2MSI/commits/master">Commits
</a></li>
<li class="hidden">
<a title="Issue Boards" class="shortcuts-issue-boards" href="/hollstei/S2MSI/board">Issue Boards</a>
</li>
</ul>
</div>

</div>
</div>
<div class="content-wrapper page-with-layout-nav">


<div class="flash-container flash-container-page">
</div>


<div class=" ">
<div class="content">
<div class="scrolling-tabs-container sub-nav-scroll">
<div class="fade-left">
<i class="fa fa-angle-left"></i>
</div>
<div class="fade-right">
<i class="fa fa-angle-right"></i>
</div>

<div class="nav-links sub-nav scrolling-tabs">
<ul class="container-fluid">
<li class="active"><a href="/hollstei/S2MSI/tree/master">Files
</a></li><li class=""><a href="/hollstei/S2MSI/commits/master">Commits
</a></li><li class=""><a href="/hollstei/S2MSI/network/master">Network
</a></li><li class=""><a href="/hollstei/S2MSI/compare?from=master&amp;to=master">Compare
</a></li><li class=""><a href="/hollstei/S2MSI/branches">Branches
</a></li><li class=""><a href="/hollstei/S2MSI/tags">Tags
</a></li></ul>
</div>
</div>

<div class="container-fluid">

<div class="tree-holder" id="tree-holder">
<div class="nav-block">
<div class="tree-ref-holder">
<form class="project-refs-form" action="/hollstei/S2MSI/refs/switch" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="destination" id="destination" value="blob" />
<input type="hidden" name="path" id="path" value="README.md" />
<div class="dropdown">
<button class="dropdown-menu-toggle js-project-refs-dropdown" type="button" data-toggle="dropdown" data-selected="master" data-ref="master" data-refs-url="/hollstei/S2MSI/refs" data-field-name="ref" data-submit-form-on-click="true"><span class="dropdown-toggle-text">master</span><i class="fa fa-chevron-down"></i></button>
<div class="dropdown-menu dropdown-menu-selectable">
<div class="dropdown-title"><span>Switch branch/tag</span><button class="dropdown-title-button dropdown-menu-close" aria-label="Close" type="button"><i class="fa fa-times dropdown-menu-close-icon"></i></button></div>
<div class="dropdown-input"><input type="search" id="" class="dropdown-input-field" placeholder="Search branches and tags" autocomplete="off" /><i class="fa fa-search dropdown-input-search"></i><i role="button" class="fa fa-times dropdown-input-clear js-dropdown-input-clear"></i></div>
<div class="dropdown-content"></div>
<div class="dropdown-loading"><i class="fa fa-spinner fa-spin"></i></div>
</div>
</div>
</form>
</div>
<ul class="breadcrumb repo-breadcrumb">
<li>
<a href="/hollstei/S2MSI/tree/master">S2MSI
</a></li>
<li>
<a href="/hollstei/S2MSI/blob/master/README.md"><strong>
README.md
</strong>
</a></li>
</ul>
</div>
<ul class="blob-commit-info hidden-xs">
<li class="commit js-toggle-container" id="commit-f2b2f625">
<a href="mailto:andre.hollstein@gfz-potsdam.de"><img class="avatar has-tooltip hidden-xs s36" alt="André Hollstein&#39;s avatar" title="André Hollstein" data-container="body" src="https://secure.gravatar.com/avatar/aac1fb736480b4b297a1f33252f31bf3?s=72&amp;d=identicon" /></a>
<div class="commit-info-block">
<div class="commit-row-title">
<span class="item-title">
<a class="commit-row-message" href="/hollstei/S2MSI/commit/f2b2f6251fb4ad9c51577b70b60bbf8664e0cd1a">Update Readme</a>
<span class="commit-row-message visible-xs-inline">
&middot;
f2b2f625
</span>
<div class="visible-xs-inline">
<a class="ci-status-link ci-status-icon-skipped " title="Commit: skipped" data-toggle="tooltip" data-placement="auto left" data-container="body" href="/hollstei/S2MSI/commit/f2b2f6251fb4ad9c51577b70b60bbf8664e0cd1a/builds"><svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 14 14">
  <g fill="#5C5C5C" fill-rule="evenodd">
    <path d="M12.5,7 C12.5,3.96243388 10.0375661,1.5 7,1.5 C3.96243388,1.5 1.5,3.96243388 1.5,7 C1.5,10.0375661 3.96243388,12.5 7,12.5 C10.0375661,12.5 12.5,10.0375661 12.5,7 Z M0,7 C0,3.13400675 3.13400675,0 7,0 C10.8659932,0 14,3.13400675 14,7 C14,10.8659932 10.8659932,14 7,14 C3.13400675,14 0,10.8659932 0,7 Z"/>
    <rect width="8" height="2" x="3" y="6" transform="rotate(45 7 7)" rx=".5"/>
  </g>
</svg>
</a>
</div>
</span>
<div class="commit-actions hidden-xs">
<a class="ci-status-link ci-status-icon-skipped " title="Commit: skipped" data-toggle="tooltip" data-placement="auto left" data-container="body" href="/hollstei/S2MSI/commit/f2b2f6251fb4ad9c51577b70b60bbf8664e0cd1a/builds"><svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 14 14">
  <g fill="#5C5C5C" fill-rule="evenodd">
    <path d="M12.5,7 C12.5,3.96243388 10.0375661,1.5 7,1.5 C3.96243388,1.5 1.5,3.96243388 1.5,7 C1.5,10.0375661 3.96243388,12.5 7,12.5 C10.0375661,12.5 12.5,10.0375661 12.5,7 Z M0,7 C0,3.13400675 3.13400675,0 7,0 C10.8659932,0 14,3.13400675 14,7 C14,10.8659932 10.8659932,14 7,14 C3.13400675,14 0,10.8659932 0,7 Z"/>
    <rect width="8" height="2" x="3" y="6" transform="rotate(45 7 7)" rx=".5"/>
  </g>
</svg>
</a>
<button class="btn btn-clipboard" data-toggle="tooltip" data-placement="bottom" data-container="body" data-clipboard-text="f2b2f6251fb4ad9c51577b70b60bbf8664e0cd1a" type="button" title="Copy to Clipboard"><i class="fa fa-clipboard"></i></button>
<a class="commit-short-id btn btn-transparent" href="/hollstei/S2MSI/commit/f2b2f6251fb4ad9c51577b70b60bbf8664e0cd1a">f2b2f625</a>

</div>
</div>
<div class="commit-row-info">
<a class="commit-author-link has-tooltip" title="andre.hollstein@gfz-potsdam.de" href="mailto:andre.hollstein@gfz-potsdam.de">André Hollstein</a>
authored
<time class="js-timeago js-timeago-pending" datetime="2016-05-03T13:08:28Z" title="May 3, 2016 3:08pm" data-toggle="tooltip" data-placement="top" data-container="body">2016-05-03 15:08:28 +0200</time><script>
//<![CDATA[
$('.js-timeago-pending').removeClass('js-timeago-pending').timeago()
//]]>
</script>
</div>
</div>
</li>

</ul>
<div class="blob-content-holder" id="blob-content-holder">
<article class="file-holder">
<div class="file-title">
<i class="fa fa-file-text-o fa-fw"></i>
<strong>
README.md
</strong>
<small>
3.71 KB
</small>
<div class="file-actions hidden-xs">
<div class="btn-group tree-btn-group">
<a class="btn btn-sm" target="_blank" href="/hollstei/S2MSI/raw/master/README.md">Raw</a>
<a class="btn btn-sm" href="/hollstei/S2MSI/blame/master/README.md">Blame</a>
<a class="btn btn-sm" href="/hollstei/S2MSI/commits/master/README.md">History</a>
<a class="btn btn-sm" href="/hollstei/S2MSI/blob/e96e2a7f28d2a6286a08bb866d356151f53754bd/README.md">Permalink</a>
</div>
<div class="btn-group" role="group">
<a class="btn btn-file-option" rel="nofollow" data-method="post" href="/hollstei/S2MSI/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fhollstei%2FS2MSI%2Fedit%2Fmaster%2FREADME.md&amp;namespace_key=86">Edit</a>
<a class="btn btn-default" rel="nofollow" data-method="post" href="/hollstei/S2MSI/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.+Try+to+replace+this+file+again.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fhollstei%2FS2MSI%2Fblob%2Fmaster%2FREADME.md&amp;namespace_key=86">Replace</a>
<a class="btn btn-remove" rel="nofollow" data-method="post" href="/hollstei/S2MSI/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.+Try+to+delete+this+file+again.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fhollstei%2FS2MSI%2Fblob%2Fmaster%2FREADME.md&amp;namespace_key=86">Delete</a>
</div>

</div>
</div>
<div class="file-content wiki">
<h1>&#x000A;<a id="introduction" class="anchor" href="#introduction" aria-hidden="true"></a>Introduction</h1>&#x000A;&#x000A;<p>The S2MSI python module is intended to allow to perform basic tasks with <a href="https://sentinel.esa.int/web/sentinel/missions/sentinel-2" rel="nofollow noreferrer" target="_blank">Sentinel-2 MSI</a> images, such as:</p>&#x000A;&#x000A;<ul>&#x000A;<li>reading tiles into numpy arrays (using either glymur, gdal, or kakadu)</li>&#x000A;<li>getting spectral response functions for a Sentinel-2</li>&#x000A;<li>manage digital elevation models for Sentinel-2 tiles</li>&#x000A;<li>getting basic data about Sentinel-2 tiles </li>&#x000A;<li>general masking capabilities (e.g. clouds and cirrus)</li>&#x000A;<li>application of machine learning tools like the classical Bayesian classifier</li>&#x000A;<li>incorporating the functionality of GRASS GIS for Sentinel-2</li>&#x000A;</ul>&#x000A;&#x000A;<p>It is also work in progress since I add functionality as I need it.</p>&#x000A;&#x000A;<h1>&#x000A;<a id="install" class="anchor" href="#install" aria-hidden="true"></a>Install</h1>&#x000A;&#x000A;<p>For now, installing is a little peculiar. I try to include functionality from GRASS GIS, for which currently only <strong>python2.7</strong> bindings exists. Apart from that, all code is <strong>python3.4+</strong>. My solution is to have everything in one module (one <strong>setup.py</strong>) which should be called with a <strong>python3.4+</strong> and a <strong>python2.7</strong> interpreter. This is automated with a bash script:</p>&#x000A;&#x000A;<pre class="code highlight js-syntax-highlight plaintext"><code>bash ./setup.sh install &#x000A;</code></pre>&#x000A;&#x000A;<p>This script expects a <strong>python2.7</strong> environment which can be activated with <code>source activate ${py27env}</code> where <code>py27env</code> is an environment variable which can be set in <code>setup.sh</code> and is set to <code>py27</code> as a default. To get everything right, the usual way of executing <code>python setup.py install</code> is not recommended as long as the GRASS GIS parts are on <strong>python2.7</strong>.</p>&#x000A;&#x000A;<p>If you want to get in touch, <a href="http://www.gfz-potsdam.de/en/section/remote-sensing/staff/profil/andre-hollstein/" rel="nofollow noreferrer" target="_blank">try this</a>. </p>&#x000A;&#x000A;<h1>&#x000A;<a id="removal" class="anchor" href="#removal" aria-hidden="true"></a>Removal</h1>&#x000A;&#x000A;<p>For now, un-installing is only possible by manually removing the folders of the module. To get a list of to-be-removed folder, call:</p>&#x000A;&#x000A;<pre class="code highlight js-syntax-highlight plaintext"><code>bash ./setup.sh uninstall&#x000A;</code></pre>&#x000A;&#x000A;<h1>&#x000A;<a id="modules" class="anchor" href="#modules" aria-hidden="true"></a>Modules</h1>&#x000A;&#x000A;<p>High level documentation is mostly missing for now, but docstrings are in place.</p>&#x000A;&#x000A;<h2>&#x000A;<a id="granuleinfo-simple-dict-like-python-object-with-spatial-information-about-sentinel-2-msi-granules" class="anchor" href="#granuleinfo-simple-dict-like-python-object-with-spatial-information-about-sentinel-2-msi-granules" aria-hidden="true"></a>GranuleInfo: Simple Dict like Python Object with Spatial Information about Sentinel-2 MSI Granules</h2>&#x000A;&#x000A;<p>Data was derived from here: <a href="https://sentinel.esa.int/web/sentinel/missions/sentinel-2/data-products" rel="nofollow noreferrer" target="_blank">https://sentinel.esa.int/web/sentinel/missions/sentinel-2/data-products</a></p>&#x000A;&#x000A;<p>There was only limited testing of the data.</p>&#x000A;&#x000A;<pre class="code highlight js-syntax-highlight python"><code><span class="kn">from</span> <span class="nn">S2MSI.GranuleInfo</span> <span class="kn">import</span> <span class="n">GranuleInfo</span>&#x000A;<span class="n">Ginfo</span> <span class="o">=</span> <span class="n">GranuleInfo</span><span class="p">(</span><span class="n">version</span><span class="o">=</span><span class="s">"lite"</span><span class="p">)</span>&#x000A;<span class="k">print</span><span class="p">(</span><span class="n">Ginfo</span><span class="p">[</span><span class="s">"32UPV"</span><span class="p">])</span>&#x000A;<span class="o">&gt;&gt;&gt;</span><span class="p">{</span><span class="s">'tr'</span><span class="p">:</span> <span class="p">{</span><span class="s">'lat'</span><span class="p">:</span> <span class="mf">49.6162737214</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">11.9045727629</span><span class="p">},</span> <span class="s">'ll'</span><span class="p">:</span> <span class="p">{</span><span class="s">'lat'</span><span class="p">:</span> <span class="mf">48.6570271781</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">10.3579107698</span><span class="p">}}</span>&#x000A;<span class="n">Ginfo</span> <span class="o">=</span> <span class="n">GranuleInfo</span><span class="p">(</span><span class="n">version</span><span class="o">=</span><span class="s">"full"</span><span class="p">)</span>&#x000A;<span class="k">print</span><span class="p">(</span><span class="n">Ginfo</span><span class="p">[</span><span class="s">"32UPV"</span><span class="p">])</span>&#x000A;<span class="o">&gt;&gt;&gt;</span><span class="p">{</span><span class="s">'pos'</span><span class="p">:</span> <span class="p">{</span><span class="s">'tl'</span><span class="p">:</span> <span class="p">{</span><span class="s">'x'</span><span class="p">:</span> <span class="mf">600000.0000025682</span><span class="p">,</span> <span class="s">'lat'</span><span class="p">:</span> <span class="mf">49.644436702</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">10.3851737332</span><span class="p">,</span> <span class="s">'y'</span><span class="p">:</span> <span class="mf">5500020.000361709</span><span class="p">},</span> <span class="s">'tr'</span><span class="p">:</span> <span class="p">{</span><span class="s">'x'</span><span class="p">:</span> <span class="mf">709800.0000165974</span><span class="p">,</span> <span class="s">'lat'</span><span class="p">:</span> <span class="mf">49.6162737214</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">11.9045727629</span><span class="p">,</span> <span class="s">'y'</span><span class="p">:</span> <span class="mf">5500020.000351718</span><span class="p">},</span> <span class="s">'lr'</span><span class="p">:</span> <span class="p">{</span><span class="s">'x'</span><span class="p">:</span> <span class="mf">709800.0000132157</span><span class="p">,</span> <span class="s">'lat'</span><span class="p">:</span> <span class="mf">48.6298215752</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">11.8474784519</span><span class="p">,</span> <span class="s">'y'</span><span class="p">:</span> <span class="mf">5390220.000321694</span><span class="p">},</span> <span class="s">'ll'</span><span class="p">:</span> <span class="p">{</span><span class="s">'x'</span><span class="p">:</span> <span class="mf">599999.9999970878</span><span class="p">,</span> <span class="s">'lat'</span><span class="p">:</span> <span class="mf">48.6570271781</span><span class="p">,</span> <span class="s">'lon'</span><span class="p">:</span> <span class="mf">10.3579107698</span><span class="p">,</span> <span class="s">'y'</span><span class="p">:</span> <span class="mf">5390220.000326163</span><span class="p">}},</span> <span class="s">'name'</span><span class="p">:</span> <span class="s">'32UPV'</span><span class="p">,</span> <span class="s">'zone'</span><span class="p">:</span> <span class="mi">32</span><span class="p">,</span> <span class="s">'epsg'</span><span class="p">:</span> <span class="mi">32632</span><span class="p">}</span>&#x000A;</code></pre>&#x000A;&#x000A;<p>Data was created with this <a href="https://git.gfz-potsdam.de/hollstei/S2MSI/tree/master/S2MSI/GranuleInfo/s2_kml_to_dict.ipynb">notebook</a>.</p>&#x000A;&#x000A;<h1>&#x000A;<a id="license" class="anchor" href="#license" aria-hidden="true"></a>License</h1>&#x000A;&#x000A;<p><a href="http://creativecommons.org/licenses/by-nc-sa/4.0/" rel="nofollow noreferrer" target="_blank"><img alt="Creative Commons License" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png"></a><br><span>S2MSI</span> by <a href="https://github.com/hollstein/S2MSI" rel="nofollow noreferrer" target="_blank">S2MSI</a> is licensed under a <a href="http://creativecommons.org/licenses/by-nc-sa/4.0/" rel="nofollow noreferrer" target="_blank">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br>Based on a work at <a href="https://github.com/hollstein/S2MSI" rel="nofollow noreferrer" target="_blank">https://github.com/hollstein/S2MSI</a>.</p>
</div>

</article>
</div>

</div>
</div>

</div>
</div>
</div>
</div>


</body>
</html>

