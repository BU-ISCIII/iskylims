"""EntryAdmin for zinnia-wymeditor"""
from django.forms import Media
from django.conf.urls import url
from django.contrib import admin
#jlgarcia 19/04/2018. Change needed for django 2.0
from django.urls import reverse
#from django.core.urlresolvers import reverse
from django.utils.translation import get_language
from django.template.response import TemplateResponse
from django.contrib.staticfiles.storage import staticfiles_storage

from zinnia.models import Entry
from zinnia.admin.entry import EntryAdmin
from zinnia.settings import ENTRY_BASE_MODEL


class EntryAdminWYMEditorMixin(object):
    """
    Mixin adding WYMeditor for editing Entry.content field.
    """

    def wymeditor(self, request):
        """
        View for serving the config of WYMEditor.
        """
        return TemplateResponse(
            request, 'admin/zinnia/entry/wymeditor.js',
            {'lang': get_language().split('-')[0]},
            content_type='application/javascript')

    def get_urls(self):
        """
        Overload the admin's urls for WYMEditor.
        """
        entry_admin_urls = super(EntryAdminWYMEditorMixin, self).get_urls()
        urls = [
            url(r'^wymeditor/$',
                self.admin_site.admin_view(self.wymeditor),
                name='zinnia_entry_wymeditor'),
        ]
        return urls + entry_admin_urls

    def _media(self):
        """
        The medias needed to enhance the admin page.
        """
        def static_url(url):
            return staticfiles_storage.url('zinnia_wymeditor/%s' % url)

        media = super(EntryAdminWYMEditorMixin, self).media

        media += Media(
            js=(static_url('js/jquery.min.js'),
                static_url('js/wymeditor/jquery.wymeditor.pack.js'),
                static_url('js/wymeditor/plugins/hovertools/'
                           'jquery.wymeditor.hovertools.js'),
                reverse('admin:zinnia_entry_wymeditor')))
        return media
    media = property(_media)


class EntryAdminWYMEditor(EntryAdminWYMEditorMixin,
                          EntryAdmin):
    """
    Enrich the default EntryAdmin with WYMEditor.
    """
    pass


if ENTRY_BASE_MODEL == 'zinnia.models_bases.entry.AbstractEntry':
    admin.site.unregister(Entry)
    admin.site.register(Entry, EntryAdminWYMEditor)
